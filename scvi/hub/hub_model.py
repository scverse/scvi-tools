import importlib
import logging
import os
from typing import Optional, Type, Union

import anndata
import torch
from anndata import AnnData
from huggingface_hub import hf_hub_download

from scvi.hub.hub_metadata import HubMetadata
from scvi.model.base import BaseModelClass

logger = logging.getLogger(__name__)


class HubModel:
    """Placeholder docstring. TODO complete"""

    def __init__(self, local_dir: str, metadata: HubMetadata):
        # TODO provide an alt ctor that takes a BaseModelClass
        self._local_dir = local_dir
        self._metadata = metadata

        self._model_path = f"{self._local_dir}/model.pt"
        self._adata_path = f"{self._local_dir}/adata.h5ad"
        self._adata_full_path = f"{self._local_dir}/adata_full.h5ad"

        # lazy load - these are not loaded until accessed
        self._model = None
        self._adata = None
        self._adata_full = None

    def push_to_huggingface_hub(self, repo_name: str, repo_url: str, repo_token: str):
        """Placeholder docstring. TODO complete"""
        raise NotImplementedError()

    @classmethod
    def pull_from_huggingface_hub(cls, repo_id: str):
        """Placeholder docstring. TODO complete"""
        m_path = hf_hub_download(repo_id=repo_id, filename="model.pt")
        d_path = hf_hub_download(repo_id=repo_id, filename="adata.h5ad")
        metadata = None
        return cls(m_path, d_path, metadata)

    def __repr__(self):
        def eval_obj(obj):
            return "No" if obj is None else "Yes"

        return (
            "HubModel with:\n"
            f"local_dir: {self._local_dir}\n"
            f"model loaded? {eval_obj(self._model)}\n"
            f"adata loaded? {eval_obj(self._adata)}\n"
            f"adata_full loaded? {eval_obj(self._adata_full)}\n"
            f"metadata:\n{self.metadata.get_formatted_attrs()}"
        )

    @property
    def metadata(self) -> HubMetadata:
        """Placeholder docstring. TODO complete"""
        return self._metadata

    @property
    def model(self) -> Type[BaseModelClass]:
        """Placeholder docstring. TODO complete"""
        if self._model is None:
            self.load_model()
        return self._model

    @property
    def adata_full(self) -> Optional[AnnData]:
        """Placeholder docstring. TODO complete"""
        if self._adata_full is None:
            self.read_adata_full()
        return self._adata_full

    @property
    def adata(self) -> AnnData:
        """Placeholder docstring. TODO complete"""
        if self._adata is None:
            self.read_adata()
        return self._adata

    def load_model(
        self,
        adata: Optional[AnnData] = None,
        use_gpu: Optional[Union[str, int, bool]] = None,
    ):
        """Placeholder docstring. TODO complete"""
        logger.info("Loading model...")
        # get the class name for this model (e.g. TOTALVI)
        torch_model = torch.load(self._model_path)
        cls_name = torch_model["attr_dict"]["registry_"]["model_name"]
        python_module = importlib.import_module("scvi.model")
        model_cls = getattr(python_module, cls_name)
        self._model = model_cls.load(
            os.path.dirname(self._model_path), adata=adata, use_gpu=use_gpu
        )

    def read_adata(self):
        """Placeholder docstring. TODO complete"""
        logger.info("Reading adata...")
        self._adata = anndata.read_h5ad(self._adata_path)

    def read_adata_full(self):
        """Download the full adata, if it exists, then read it into memory."""
        # logger.info("Reading adata_full...")
        # self._adata_full = anndata.read_h5ad(self._adata_full_path)
        # -> need to download it using the metadata url (or wherever it is) first
        raise NotImplementedError()
