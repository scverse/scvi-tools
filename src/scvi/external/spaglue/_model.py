import logging
from itertools import cycle

import numpy as np
import scipy.sparse
import torch
from anndata import AnnData
from torch.utils.data import DataLoader
from torch_geometric.data import Data

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
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

        generative_distributions = {
            mod: adata.uns["spaglue_likelihood"] for mod, adata in self.adatas.items()
        }

        if guidance_graph is not None:
            self.guidance_graph = guidance_graph
        else:
            self.guidance_graph = _construct_guidance_graph(self.adatas)

        _check_guidance_graph_consisteny(self.guidance_graph, adatas)

        self.module = SPAGLUEVAE(
            n_inputs=n_inputs,
            n_batches=n_batches,
            gene_likelihoods=generative_distributions,
            guidance_graph=self.guidance_graph,
            # **model_kwargs,
        )

        self._model_summary_string = (
            f"SpaGlue Model with the following params: modalities: {self.modality_names}, "
            f"n_inputs: {n_inputs}, n_batches: {n_batches}, "
            f"generative distributions: {generative_distributions}"
        )

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
        # **kwargs: dict,
    ) -> None:
        if scipy.sparse.issparse(adata.X) and not isinstance(adata.X, scipy.sparse.csr_matrix):
            adata.X = adata.X.tocsr()

        # For a specific layer (e.g., "counts")
        if "counts" in adata.layers and not isinstance(
            adata.layers["counts"], scipy.sparse.csr_matrix
        ):
            adata.layers["counts"] = adata.layers["counts"].tocsr()

        adata.uns["spaglue_likelihood"] = likelihood

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
