from collections.abc import Iterator, Sequence
from itertools import cycle

import numpy as np
import torch
from anndata import AnnData
from torch.utils.data import DataLoader

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalObsField, LayerField
from scvi.dataloaders import DataSplitter
from scvi.model._utils import parse_device_args
from scvi.model.base import BaseModelClass
from scvi.train import Trainer

from ._module import SPAGLUEVAE
from ._task import SPAGLUETrainingPlan


class SPAGLUE(BaseModelClass):
    def __init__(
        self,
        adata_seq: AnnData,
        adata_spatial: AnnData,
        generative_distributions: list[str] | None = None,
        # **model_kwargs: dict,
    ) -> None:
        super().__init__()
        self.adatas = [adata_seq, adata_spatial]

        self.adata_managers = {
            "seq": self._get_most_recent_anndata_manager(adata_seq, required=True),
            "spatial": self._get_most_recent_anndata_manager(adata_spatial, required=True),
        }

        self.registries_ = []
        # register data_manager for this specific instance
        for adm in self.adata_managers.values():
            self._register_manager_for_instance(adm)
            self.registries_.append(adm.registry)

        sum_stats = [adm.summary_stats for adm in self.adata_managers.values()]

        n_inputs = [s["n_vars"] for s in sum_stats]
        n_batches = [s["n_batch"] for s in sum_stats]

        generative_distributions = generative_distributions or ["nb", "nb"]  ## for now

        self.module = SPAGLUEVAE(
            n_inputs=n_inputs,
            n_batches=n_batches,
            gene_likelihoods=generative_distributions,
            # **model_kwargs,
        )

        self._model_summary_string = (
            f"SpaGlue Model with the following params: n_inputs_seq: {n_inputs[0]}, "
            f"n_inputs_spa: {n_inputs[1]} , n_batches_seq: {n_batches[0]}, "
            f"n_batches_spa: {n_batches[1]} , generative distributions: {generative_distributions}"
        )
        self._module_class = SPAGLUEVAE

    def train(
        self,
        max_epochs: int = 200,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float = 0.9,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        # datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        # **kwargs: dict,
    ) -> None:
        if torch.backends.mps.is_available():
            accelerator = "mps"
            devices = 1  # Use a single MPS device
            print("Using MPS accelerator for training.")
        else:
            print("MPS not available. Falling back to default accelerator.")

        accelerator, devices, device = parse_device_args(
            accelerator=accelerator,  # cpu, gpu or auto (automatically selects optimal device)
            devices=devices,  # auto, 1 0r [0,1]
            return_device="torch",  #  make returned device pytorch compatible
        )

        # datasplitter_kwargs = datasplitter_kwargs or {}

        self.trainer = Trainer(
            max_epochs=max_epochs,
            enable_progress_bar=True,
            log_every_n_steps=1,
            early_stopping=True,
            early_stopping_monitor="validation_loss",
            accelerator=accelerator,
            devices=devices,
            # **kwargs,
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

        train_dl = TrainDL(train_dls)  # combine the list of TRAINING dataloaders
        val_dl = TrainDL(val_dls)

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else {}
        self._training_plan = SPAGLUETrainingPlan(
            self.module,  # JVAE module (defined in __init__)
            **plan_kwargs,
        )

        if train_size == 1.0:  # avoid issue with empty train/val dataloaders
            # circumvent the empty data loader problem if all dataset used for training
            self.trainer.fit(
                self._training_plan, train_dl
            )  # calling trainer.fit uses train_df and val_dl
        else:
            # accepts list of val dataloaders
            self.trainer.fit(
                self._training_plan, train_dataloaders=train_dl, val_dataloaders=val_dl
            )

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

    @torch.inference_mode()
    def get_latent_representation(
        self,
        adata_seq: AnnData | None = None,
        adata_spatial: AnnData | None = None,
        indices_seq: Sequence[int] | None = None,
        indices_spatial: Sequence[int] | None = None,
        batch_size: int | None = None,
    ) -> dict[str, np.ndarray]:
        self._check_if_trained(warn=False)

        results = {}

        # Compute latent space for sequencing data if provided
        if adata_seq is not None:
            adata_seq = self._validate_anndata(adata_seq)
            dataloader_seq = self._make_data_loader(
                adata=adata_seq, indices=indices_seq, batch_size=batch_size
            )
            results["seq"] = self._compute_latent(dataloader_seq, mode=0)

        # Compute latent space for spatial data if provided
        if adata_spatial is not None:
            adata_spatial = self._validate_anndata(adata_spatial)
            dataloader_spatial = self._make_data_loader(
                adata=adata_spatial, indices=indices_spatial, batch_size=batch_size
            )
            results["spatial"] = self._compute_latent(dataloader_spatial, mode=1)

        return results

    def _compute_latent(
        self,
        dataloader: Iterator[dict[str, torch.Tensor | None]],
        mode: int,
    ) -> np.ndarray:
        zs: list[torch.Tensor] = []

        for batch in dataloader:
            latent_tensor = self.module.inference(
                **self.module._get_inference_input(batch), mode=mode
            )["z"]

            latent_tensor = latent_tensor.cpu().numpy()

            zs.append(latent_tensor)

        return np.concatenate(zs, axis=0)


class TrainDL(DataLoader):  # creates batch structure for training process
    # def __init__(self, data_loader_list, **kwargs):
    def __init__(self, data_loader_list):
        self.data_loader_list = data_loader_list  # list of individual dls
        self.largest_train_dl_idx = np.argmax(
            [len(dl.indices) for dl in data_loader_list]
        )  # index of dl with largest num of samples
        self.largest_dl = self.data_loader_list[
            self.largest_train_dl_idx
        ]  # dl corresponding to the largest dataset
        # super().__init__(self.largest_dl, **kwargs)
        super().__init__(self.largest_dl)

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
