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
from scvi.model._utils import parse_device_args
from scvi.model.base import BaseModelClass, VAEMixin
from scvi.module._constants import MODULE_KEYS
from scvi.train import Trainer

from ._module import SPAGLUEVAE
from ._task import SPAGLUETrainingPlan

logger = logging.getLogger(__name__)


class SPAGLUE(BaseModelClass, VAEMixin):
    def __init__(
        self,
        adata_seq: AnnData,
        adata_spatial: AnnData,
        generative_distributions: list[str] | None = None,
        guidance_graph: Data | None = None,
        **model_kwargs: dict,
    ) -> None:
        super().__init__()
        self.adatas = [adata_seq, adata_spatial]

        self.adata_managers = {
            "seq": self._get_most_recent_anndata_manager(adata_seq, required=True),
            "spatial": self._get_most_recent_anndata_manager(adata_spatial, required=True),
        }

        self.registries_ = []
        for adm in self.adata_managers.values():
            self._register_manager_for_instance(adm)
            self.registries_.append(adm.registry)

        sum_stats = [adm.summary_stats for adm in self.adata_managers.values()]

        n_inputs = [s["n_vars"] for s in sum_stats]
        n_batches = [s["n_batch"] for s in sum_stats]

        generative_distributions = generative_distributions or ["nb", "nb"]  ## for now

        if guidance_graph is not None:
            self.guidance_graph = guidance_graph
        else:
            self.guidance_graph = _construct_guidance_graph(adata_seq, adata_spatial)

        _check_guidance_graph_consisteny(self.guidance_graph, self.adatas)

        self.module = SPAGLUEVAE(
            n_inputs=n_inputs,
            n_batches=n_batches,
            gene_likelihoods=generative_distributions,
            guidance_graph=self.guidance_graph,
            # **model_kwargs,
        )

        self._model_summary_string = (
            f"SpaGlue Model with the following params: n_inputs_seq: {n_inputs[0]}, "
            f"n_inputs_spa: {n_inputs[1]} , n_batches_seq: {n_batches[0]}, "
            f"n_batches_spa: {n_batches[1]} , generative distributions: {generative_distributions}"
        )
        # self._module_class = SPAGLUEVAE

    def train(
        self,
        max_epochs: int = 200,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float = 0.9,
        shuffle_set_split: bool = True,
        batch_size: int = 256,
        # datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,  # kwargs passed to trainingplan
        # graph_loss_weight = self.graph_loss_weight,
        **kwargs,
    ) -> None:
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
        train_dls, test_dls, val_dls = [], [], []
        for _i, adm in enumerate(
            self.adata_managers.values()
        ):  # one adata manager for each modality
            ds = DataSplitter(  # creates dataloaders as part of ds setup process
                adm,
                train_size=train_size,  # train_size = None defaults to 0.9
                validation_size=validation_size,
                batch_size=batch_size,
                shuffle_set_split=shuffle_set_split,
                # **datasplitter_kwargs,
            )
            ds.setup()
            train_dls.append(ds.train_dataloader())
            test_dls.append(ds.test_dataloader())
            val_dls.append(ds.val_dataloader())

            self.train_indices_.append(ds.train_idx)
            self.test_indices_.append(ds.test_idx)
            self.validation_indices_.append(ds.val_idx)

        # train_dl = TrainDL(train_dls, num_workers=9)  # combine the list of TRAINING dataloaders
        # val_dl = TrainDL(val_dls, num_workers=9)

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
        self.module.eval()

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
        # **kwargs: dict,
    ) -> None:
        if scipy.sparse.issparse(adata.X) and not isinstance(adata.X, scipy.sparse.csr_matrix):
            adata.X = adata.X.tocsr()

        # For a specific layer (e.g., "counts")
        if "counts" in adata.layers and not isinstance(
            adata.layers["counts"], scipy.sparse.csr_matrix
        ):
            adata.layers["counts"] = adata.layers["counts"].tocsr()

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
        self, adatas: list[AnnData] = None, batch_size: int = 128
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
        adatas: list[AnnData] = None,
        deterministic: bool = True,
        batch_size: int = 128,
    ) -> list[np.ndarray]:
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
        scdls = self._make_scvi_dls(adatas, batch_size=batch_size)
        self.module.eval()
        latents = []
        for mode, scdl in enumerate(scdls):
            latent = []
            for tensors in scdl:
                latent.append(
                    self.module.inference(**self.module._get_inference_input(tensors), mode=mode)[
                        "z"
                    ]
                    .cpu()
                    .detach()
                )

            latent = torch.cat(latent).numpy()
            latents.append(latent)

        return latents

    @torch.inference_mode()
    def get_imputed_values(
        self,
        source_mode: int,
        source_adata: AnnData | None = None,
        batch_size: int = 128,
        target_batch: int | None = None,
        target_libsize: float | None = None,
    ) -> list[np.ndarray]:
        # choose source adata according to mode
        if source_mode not in (0, 1):
            raise ValueError("`source_mode` must be 0 or 1!")

        if source_adata is None:
            source_adata = self.adatas[source_mode]

        self.module.eval()
        dl = self._make_data_loader(source_adata, batch_size=batch_size)

        reconstructed_counts = []
        for tensor in dl:
            inference_output = self.module.inference(
                **self.module._get_inference_input(tensor), mode=source_mode
            )

            # keep all spatial except the embedded features
            inference_output["v"] = inference_output["v_other"]

            generative_input = self.module._get_generative_input(
                tensor,
                inference_output,
            )

            # choose batch 0 by default
            if target_batch is None:
                target_batch = 0

            if not isinstance(target_batch, int):
                raise TypeError("`target_batch` must be an integer.")

            generative_input[MODULE_KEYS.BATCH_INDEX_KEY] = torch.full(
                (generative_input[MODULE_KEYS.LIBRARY_KEY].shape[0], 1),  # batch_size x 1
                fill_value=target_batch,
                dtype=torch.long,
                device=generative_input[MODULE_KEYS.LIBRARY_KEY].device,
            )

            # Overwrite library size if needed
            if target_libsize is not None:
                if not isinstance(target_libsize, float):
                    raise TypeError("`target_libsize` must be an integer.")

                generative_input[MODULE_KEYS.LIBRARY_KEY] = torch.tensor(
                    target_libsize,
                    dtype=generative_input[MODULE_KEYS.LIBRARY_KEY].dtype,
                    device=generative_input[MODULE_KEYS.LIBRARY_KEY].device,
                )

            target_mode = 1 - source_mode

            generative_output = self.module.generative(**generative_input, mode=target_mode)

            reconstructed_counts.append(generative_output["px_rate"].cpu().detach())

        reconstructed_count = torch.cat(reconstructed_counts).numpy()
        return reconstructed_count


class TrainDL(DataLoader):  # creates batch structure for training process
    def __init__(self, data_loader_list):
        # def __init__(self, data_loader_list, num_workers=0):
        self.data_loader_list = data_loader_list  # list of individual dls
        self.largest_train_dl_idx = np.argmax(
            [len(dl.indices) for dl in data_loader_list]
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
        return zip(
            *train_dls, strict=False
        )  # zips iterators together - one batch for each dataset each time
        # until larger one runs out
