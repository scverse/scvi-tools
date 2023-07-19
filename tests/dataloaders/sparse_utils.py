from typing import List, Literal, Optional, Union

import torch
from anndata import AnnData

import scvi
from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager, fields
from scvi.train import TrainRunner


class TestSparseDataSplitter(scvi.dataloaders.DataSplitter):
    def __init__(
        self, *args, expected_sparse_layout: Literal["csr", "csc"] = None, **kwargs
    ):
        if expected_sparse_layout == "csr":
            self.expected_sparse_layout = torch.sparse_csr
        elif expected_sparse_layout == "csc":
            self.expected_sparse_layout = torch.sparse_csc
        else:
            self.expected_sparse_layout = None

        super().__init__(*args, **kwargs)

    def on_after_batch_transfer(self, batch, dataloader_idx):
        X = batch.get(scvi.REGISTRY_KEYS.X_KEY)
        assert isinstance(X, torch.Tensor)
        assert X.layout is self.expected_sparse_layout

        batch = super().on_after_batch_transfer(batch, dataloader_idx)

        X = batch.get(scvi.REGISTRY_KEYS.X_KEY)
        assert isinstance(X, torch.Tensor)
        assert X.layout == torch.strided

        return batch


class TestSparseTrainingPlan(scvi.train.TrainingPlan):
    def training_step(self, batch, batch_idx):
        pass

    def validation_step(self, batch, batch_idx):
        pass


class TestSparseModule(scvi.module.base.BaseModuleClass):
    def __init__(self):
        super().__init__()
        self.layer = torch.nn.Linear(1, 1)

    def _get_inference_input(self):
        pass

    def _get_generative_input(self):
        pass

    def inference(self):
        pass

    def generative(self):
        pass

    def loss(self):
        pass


class TestSparseModel(scvi.model.base.BaseModelClass):
    def __init__(self, adata: AnnData):
        super().__init__(adata)
        self.module = TestSparseModule()

    @classmethod
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: Optional[str] = None,
        batch_key: Optional[str] = None,
    ):
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            fields.LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            fields.CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata)
        cls.register_manager(adata_manager)

    def train(
        self,
        max_epochs: int = 1,
        accelerator: str = "auto",
        devices: Union[int, List[int], str] = "auto",
        expected_sparse_layout: Literal["csr", "csc"] = None,
    ):
        data_splitter = TestSparseDataSplitter(
            self.adata_manager,
            expected_sparse_layout=expected_sparse_layout,
            load_sparse_tensor=True,
        )
        training_plan = TestSparseTrainingPlan(self.module)
        runner = TrainRunner(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
        )
        return runner()
