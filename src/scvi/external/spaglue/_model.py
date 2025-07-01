import logging
import os
import warnings
from itertools import cycle
from typing import Literal

import anndata as ad
import numpy as np
import scipy.sparse
import torch
from anndata import AnnData
from torch.utils.data import DataLoader
from torch_geometric.data import Data

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._constants import _MODEL_NAME_KEY, _SETUP_ARGS_KEY
from scvi.data._download import _download
from scvi.data.fields import CategoricalObsField, LayerField
from scvi.dataloaders import AnnDataLoader, DataSplitter
from scvi.external.spaglue._utils import (
    _check_guidance_graph_consisteny,
    _construct_guidance_graph,
)
from scvi.model._utils import get_max_epochs_heuristic, parse_device_args
from scvi.model.base import BaseModelClass, VAEMixin
from scvi.module._constants import MODULE_KEYS
from scvi.train import Trainer

from ._module import SPAGLUEVAE
from ._task import SPAGLUETrainingPlan

logger = logging.getLogger(__name__)


def _load_saved_spaglue_files(
    dir_path: str,
    prefix: str | None = None,
    map_location: Literal["cpu", "cuda"] | None = None,
    backup_url: str | None = None,
) -> tuple[dict, dict, np.ndarray, np.ndarray, dict, AnnData | None, AnnData | None]:
    file_name_prefix = prefix or ""

    model_file_name = f"{file_name_prefix}model.pt"
    model_path = os.path.join(dir_path, model_file_name)

    try:
        _download(backup_url, dir_path, model_file_name)
        model = torch.load(model_path, map_location=map_location, weights_only=False)
    except FileNotFoundError as exc:
        raise ValueError(f"Failed to load model file at {model_path}. ") from exc

    modality_names = model["modality_names"]

    adatas = {}
    var_names = {}
    for mod in modality_names:
        adata_path = os.path.join(dir_path, f"{file_name_prefix}adata_{mod}.h5ad")
        if os.path.exists(adata_path):
            adatas[mod] = ad.read_h5ad(adata_path)
            var_names[mod] = adatas[mod].var_names
        else:
            adatas[mod] = None

    model_state_dict = model["model_state_dict"]
    attr_dict = model["attr_dict"]

    return (
        attr_dict,
        var_names,
        model_state_dict,
        adatas,
    )


class SPAGLUE(BaseModelClass, VAEMixin):
    def __init__(
        self,
        adatas: dict[str, AnnData],
        guidance_graph: Data | None = None,
        **model_kwargs: dict,
    ) -> None:
        super().__init__()
        self.adatas = adatas
        self.modality_names = list(adatas.keys())

        self.adata_managers = {
            mod: self._get_most_recent_anndata_manager(adata, required=True)
            for mod, adata in self.adatas.items()
        }

        self.registries_ = []
        for adm in self.adata_managers.values():
            self._register_manager_for_instance(adm)
            self.registries_.append(adm.registry)

        sum_stats = {mod: adm.summary_stats for mod, adm in self.adata_managers.items()}
        self.summary_stats = sum_stats

        n_inputs = {mod: s["n_vars"] for mod, s in sum_stats.items()}
        n_batches = {mod: s["n_batch"] for mod, s in sum_stats.items()}
        n_labels = {mod: s["n_labels"] for mod, s in sum_stats.items()}

        # Which generative distributions to use
        generative_distributions = {
            mod: adata.uns["spaglue_likelihood"] for mod, adata in self.adatas.items()
        }

        # Whether to use the Gaussian Mixture Prior
        gmm_priors = {mod: adata.uns["spaglue_gmm_prior"] for mod, adata in self.adatas.items()}

        # Whether to use the provided label key to guide model development
        semi_supervised = {
            mod: adata.uns["spaglue_semi_supervised"] for mod, adata in self.adatas.items()
        }

        # How many components to model
        n_mixture_components = {
            mod: adata.uns["spaglue_n_mixture_components"] for mod, adata in self.adatas.items()
        }

        # Compute guidance graph if not provided, do a sanity check
        if guidance_graph is not None:
            self.guidance_graph = guidance_graph
        else:
            self.guidance_graph = _construct_guidance_graph(self.adatas)
        _check_guidance_graph_consisteny(self.guidance_graph, adatas)

        self.module = SPAGLUEVAE(
            n_inputs=n_inputs,
            n_batches=n_batches,
            n_labels=n_labels,
            gene_likelihoods=generative_distributions,
            guidance_graph=self.guidance_graph,
            use_gmm_prior=gmm_priors,
            semi_supervised=semi_supervised,
            n_mixture_components=n_mixture_components,
            **model_kwargs,
        )

        self._model_summary_string = (
            f"SpaGlue Model with the following params: modalities: {self.modality_names}, "
            f"n_inputs: {n_inputs}, n_batches: {n_batches}, "
            f"generative distributions: {generative_distributions}"
        )
        self.init_params_ = self._get_init_params(locals())

    def train(
        self,
        max_epochs: int = 200,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float = 0.9,
        shuffle_set_split: bool = True,
        batch_size: int = 1024,
        # datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,  # kwargs passed to trainingplan
        **kwargs,
    ) -> None:
        if max_epochs is None:
            min_obs = np.min(
                [self.summary_stats["diss"]["n_cells"], self.summary_stats["spatial"]["n_cells"]]
            )
            max_epochs = get_max_epochs_heuristic(min_obs)
            logger.info(f"max_epochs was approximated to {max_epochs}")

        accelerator, devices, device = parse_device_args(
            accelerator=accelerator,  # cpu, gpu or auto (automatically selects optimal device)
            devices=devices,  # auto, 1 0r [0,1]
            return_device="torch",  #  make returned device pytorch compatible
        )

        self.trainer = Trainer(
            max_epochs=max_epochs,
            enable_progress_bar=True,
            progress_bar_refresh_rate=1,
            early_stopping=True,
            early_stopping_monitor="validation_loss",
            early_stopping_min_delta=1,
            early_stopping_patience=10,
            early_stopping_mode="min",
            accelerator=accelerator,
            devices=devices,
            **kwargs,
        )

        validation_size = 1 - train_size
        self.train_indices_, self.test_indices_, self.validation_indices_ = [], [], []
        train_dls, test_dls, val_dls = {}, {}, {}
        for mod, adm in self.adata_managers.items():  # mod is the modality name
            ds = DataSplitter(
                adm,
                train_size=train_size,
                validation_size=validation_size,
                batch_size=batch_size,
                shuffle_set_split=shuffle_set_split,
            )
            ds.setup()
            train_dls[mod] = ds.train_dataloader()
            test_dls[mod] = ds.test_dataloader()
            val_dls[mod] = ds.val_dataloader()

            self.train_indices_.append(ds.train_idx)
            self.test_indices_.append(ds.test_idx)
            self.validation_indices_.append(ds.val_idx)

        train_dl = TrainDL(train_dls)  # combine the list of TRAINING dataloaders
        val_dl = TrainDL(val_dls)

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else {}

        self._training_plan = SPAGLUETrainingPlan(
            self.module,  # JVAE module (defined in __init__)
            **plan_kwargs,
        )

        if train_size == 1.0:  # avoid issue with empty train/val dataloaders
            self.trainer.fit(
                self._training_plan, train_dl
            )  # calling trainer.fit uses train_df and val_dl
        else:
            self.trainer.fit(
                self._training_plan, train_dataloaders=train_dl, val_dataloaders=val_dl
            )
        try:
            self.history_ = self.trainer.logger.history
        except AttributeError:
            self.history_ = None

        self.module.eval()  # set model to evaluation mode (di)

        self.to_device(device)  # move model to speciified device
        self.is_trained_ = True  # marl model as trained

    @classmethod
    def setup_anndata(
        cls,
        adata: AnnData,
        batch_key: str | None = None,
        labels_key: str | None = None,
        layer: str | None = None,
        likelihood: str = "nb",
        gmm_prior: bool = False,
        semi_supervised: bool = False,
        n_mixture_components: int = 10,
        **kwargs: dict,
    ) -> None:
        if scipy.sparse.issparse(adata.X) and not isinstance(adata.X, scipy.sparse.csr_matrix):
            adata.X = adata.X.tocsr()

        # For a specific layer (e.g., "counts")
        if "counts" in adata.layers and not isinstance(
            adata.layers["counts"], scipy.sparse.csr_matrix
        ):
            adata.layers["counts"] = adata.layers["counts"].tocsr()

        adata.uns["spaglue_likelihood"] = likelihood
        adata.uns["spaglue_gmm_prior"] = gmm_prior
        adata.uns["spaglue_semi_supervised"] = semi_supervised
        if semi_supervised:
            adata.uns["spaglue_n_mixture_components"] = len(adata.obs[labels_key].unique())
        else:
            adata.uns["spaglue_n_mixture_components"] = n_mixture_components

        # Set up the anndata object for the model
        setup_method_args = cls._get_setup_method_args(
            **locals()
        )  # returns dict organizing the args used to call setup anndata
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        # adata_manager.register_fields(adata, **kwargs)
        adata_manager.register_fields(adata)
        cls.register_manager(adata_manager)

    def _make_scvi_dls(
        self, adatas: list[AnnData] = None, batch_size: int = 1024
    ) -> list[AnnDataLoader]:
        if adatas is None:
            adatas = self.adatas
        post_list = [self._make_data_loader(ad, batch_size=batch_size) for ad in adatas]
        for i, dl in enumerate(post_list):
            dl.mode = i

        return post_list

    @torch.inference_mode()
    def get_latent_representation(
        self,
        # adatas: list[AnnData] = None,
        adatas: dict[str, AnnData] | list[AnnData] = None,
        deterministic: bool = True,
        batch_size: int = 1024,
    ) -> dict[np.ndarray]:
        """Return the latent space embedding for each dataset.

        Parameters
        ----------
        adatas
            List of adata seq and adata spatial.
        deterministic
            If true, use the mean of the encoder instead of a Gaussian sample.
        batch_size
            Minibatch size for data loading into model.
        """
        if adatas is None:
            adatas = self.adatas
        if isinstance(adatas, dict):
            modality_names = list(adatas.keys())
            adata_list = list(adatas.values())
        else:
            modality_names = self.modality_names
            adata_list = adatas

        scdls = self._make_scvi_dls(adata_list, batch_size=batch_size)

        self.module.eval()
        latents = {}
        for modality, scdl in zip(modality_names, scdls, strict=False):
            latent = []
            for tensors in scdl:
                latent.append(
                    self.module.inference(
                        **self.module._get_inference_input(tensors), mode=modality
                    )["z"]
                    .cpu()
                    .detach()
                )
            latent = torch.cat(latent).numpy()
            latents[modality] = latent

        return latents

    @torch.inference_mode()
    def get_imputed_values(
        self,
        source_modality: int,
        source_adata: AnnData | None = None,
        batch_size: int = 1024,
        target_batch: int | None = None,
        target_libsize: float | None = None,
    ) -> list[np.ndarray]:
        # choose source adata according to mode
        if source_modality not in self.adatas:
            raise ValueError(f"`source_modality` must be one of {list(self.adatas.keys())}!")

        # Choose source AnnData
        if source_adata is None:
            source_adata = self.adatas[source_modality]

        # Get batch categories for this modality
        batch_manager = self.adata_managers[source_modality]
        batch_categories = batch_manager.get_state_registry(REGISTRY_KEYS.BATCH_KEY).get(
            "categorical_mapping", None
        )

        self.module.eval()
        dl = self._make_data_loader(source_adata, batch_size=batch_size)

        reconstructed_counts = []
        for tensor in dl:
            inference_output = self.module.inference(
                **self.module._get_inference_input(tensor), mode=source_modality
            )

            # keep all spatial except the embedded features
            inference_output["v"] = inference_output["v_other"]
            generative_input = self.module._get_generative_input(
                tensor,
                inference_output,
            )

            # --- Handle target batch ---
            batch_size_ = generative_input[MODULE_KEYS.LIBRARY_KEY].shape[0]
            if target_batch is not None:
                b = np.asarray(target_batch)
                if b.size == 1:  # if size 1: broadcast to array of length Batch size
                    b = np.full(batch_size_, b.item())
                elif b.size != batch_size_:  # raise error if wrong size
                    raise ValueError("`target_batch` must have the same size as adata!")
                # Map categorical names to indices if needed
                if batch_categories is not None and not np.issubdtype(b.dtype, np.integer):
                    b = np.array([np.where(batch_categories == lbl)[0][0] for lbl in b])
            # use batch index zero if no target batch is provided
            else:
                b = np.zeros(batch_size_, dtype=int)
            b = np.asarray(b, dtype=int).reshape(-1, 1)

            generative_input[MODULE_KEYS.BATCH_INDEX_KEY] = torch.tensor(
                b,
                dtype=torch.long,
                device=generative_input[MODULE_KEYS.LIBRARY_KEY].device,
            )

            # --- Handle target libsize ---
            if target_libsize is not None:
                l = target_libsize
                if not isinstance(l, np.ndarray):
                    l = np.asarray(l)
                l = l.squeeze()

                if l.ndim == 0:  # Scalar
                    l = l[np.newaxis]
                elif l.ndim > 1:
                    raise ValueError("`target_libsize` cannot be >1 dimensional")

                if l.size == 1:
                    l = np.repeat(l, batch_size_)

                if l.size != batch_size_:
                    raise ValueError("`target_libsize` must have the same size as the batch!")

                l = l.reshape((-1, 1))
                generative_input[MODULE_KEYS.LIBRARY_KEY] = torch.tensor(
                    l,
                    dtype=generative_input[MODULE_KEYS.LIBRARY_KEY].dtype,
                    device=generative_input[MODULE_KEYS.LIBRARY_KEY].device,
                )

            # Determine target modality (the other one)
            target_modalities = [m for m in self.modality_names if m != source_modality]
            if len(target_modalities) != 1:
                raise ValueError("There must be exactly two modalities defined.")
            target_modality = target_modalities[0]

            generative_output = self.module.generative(**generative_input, mode=target_modality)

            reconstructed_counts.append(generative_output["px_rate"].cpu().detach())

        reconstructed_count = torch.cat(reconstructed_counts).numpy()
        return reconstructed_count

    def save(
        self,
        dir_path: str,
        prefix: str | None = None,
        overwrite: bool = False,
        save_anndata: bool = False,
        save_kwargs: dict | None = None,
        **anndata_write_kwargs,
    ):
        if not os.path.exists(dir_path) or overwrite:
            os.makedirs(dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                f"{dir_path} already exists. Please provide an unexisting directory for saving."
            )

        file_name_prefix = prefix or ""
        save_kwargs = save_kwargs or {}

        adatas = self.adatas
        if save_anndata:
            for key in self.modality_names:
                ad = adatas[key]
                save_path = os.path.join(dir_path, f"adata_{key}.h5ad")
                ad.write(save_path)

        model_state_dict = self.module.state_dict()

        var_names = {}
        for key in self.modality_names:
            var_names[key] = adatas[key].var_names

        # get all the user attributes
        user_attributes = self._get_user_attributes()

        # only save the public attributes with _ at the very end
        user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}

        model_save_path = os.path.join(dir_path, f"{file_name_prefix}model.pt")

        torch.save(
            {
                "model_state_dict": model_state_dict,
                "modality_names": self.modality_names,
                **{f"var_names_{mod}": var_names[mod] for mod in self.modality_names},
                "attr_dict": user_attributes,
            },
            model_save_path,
            **save_kwargs,
        )

    @classmethod
    # @devices_dsp.dedent
    def load(
        cls,
        dir_path: str,
        adata_seq: AnnData | None = None,
        adata_spatial: AnnData | None = None,
        accelerator: str = "auto",
        device: int | str = "auto",
        prefix: str | None = None,
        backup_url: str | None = None,
    ):
        _, _, device = parse_device_args(
            accelerator=accelerator,
            devices=device,
            return_device="torch",
        )

        (
            attr_dict,
            var_names_dict,
            model_state_dict,
            adatas,
        ) = _load_saved_spaglue_files(
            dir_path,
            prefix=prefix,
            map_location=device,
            backup_url=backup_url,
        )

        for mod in adatas.keys():
            saved_var_names = var_names_dict[mod]
            user_var_names = adatas[mod].var_names.astype(str)

            if not np.array_equal(saved_var_names, user_var_names):
                warnings.warn(
                    "var_names for adata passed in does not match var_names of adata "
                    "used to train the model. For valid results, the vars need to be the"
                    "same and in the same order as the adata used to train the model.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )

        registries = attr_dict.pop("registries_")
        for adata, registry in zip(adatas, registries, strict=True):
            if _MODEL_NAME_KEY in registry and registry[_MODEL_NAME_KEY] != cls.__name__:
                raise ValueError("It appears you are loading a model from a different class.")

            if _SETUP_ARGS_KEY not in registry:
                raise ValueError(
                    "Saved model does not contain original setup inputs. "
                    "Cannot load the original setup."
                )

            cls.setup_anndata(adatas[adata], source_registry=registry, **registry[_SETUP_ARGS_KEY])

        # get the parameters for the class init signature
        init_params = attr_dict.pop("init_params_")

        # new saving and loading, enable backwards compatibility
        if "non_kwargs" in init_params.keys():
            # grab all the parameters except for kwargs (is a dict)
            non_kwargs = init_params["non_kwargs"]
            kwargs = init_params["kwargs"]

            # expand out kwargs
            kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        else:
            # grab all the parameters except for kwargs (is a dict)
            non_kwargs = {k: v for k, v in init_params.items() if not isinstance(v, dict)}
            kwargs = {k: v for k, v in init_params.items() if isinstance(v, dict)}
            kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}

        # Remove 'adatas' from non_kwargs and kwargs if present
        non_kwargs.pop("adatas", None)
        kwargs.pop("adatas", None)

        model = cls(adatas, **non_kwargs, **kwargs)

        for attr, val in attr_dict.items():
            setattr(model, attr, val)

        model.module.load_state_dict(model_state_dict)
        model.module.eval()
        model.to_device(device)
        return model


class TrainDL(DataLoader):  # creates batch structure for training process
    def __init__(self, data_loader_dict):
        # def __init__(self, data_loader_list, num_workers=0):
        self.data_loader_dict = data_loader_dict
        self.modality_names = list(data_loader_dict.keys())
        self.data_loader_list = list(data_loader_dict.values())

        self.largest_train_dl_idx = np.argmax(
            [len(dl.indices) for dl in self.data_loader_list]
        )  # index of dl with largest num of samples
        self.largest_dl = self.data_loader_list[
            self.largest_train_dl_idx
        ]  # dl corresponding to the largest dataset
        # super().__init__(self.largest_dl, num_workers=num_workers)
        super().__init__(
            self.largest_dl,
            num_workers=settings.dl_num_workers,
            persistent_workers=getattr(settings, "dl_persistent_workers", False),
        )

    # number of batches per epoch is determined by the larger dataset
    def __len__(self):
        return len(self.largest_dl)

    # cyclic iteration
    def __iter__(self):
        # produces batches by zipping together
        # ensures that smaller datasets are repeated indefinitely
        train_dls = [
            (
                dl if i == self.largest_train_dl_idx else cycle(dl)
            )  # repeat smaller dataset indefinitely
            for i, dl in enumerate(self.data_loader_list)
        ]
        # Return a dict of modality_name -> batch for each step
        for batches in zip(*train_dls, strict=False):
            yield dict(zip(self.modality_names, batches, strict=False))
