import importlib
import logging
import os
from pathlib import Path
from typing import List, Optional, Type, Union

import anndata
import torch
from anndata import AnnData
from huggingface_hub import (
    HfApi,
    ModelCard,
    ModelFilter,
    create_repo,
    snapshot_download,
)
from huggingface_hub.hf_api import ModelInfo

from scvi.model.base import BaseModelClass

HF_LIBRARY_NAME = "scvi-tools"
logger = logging.getLogger(__name__)


class HubModel:
    """Placeholder docstring. TODO complete"""

    def __init__(
        self,
        local_dir: str,
        model_card: Union[ModelCard, str],
    ):
        self._local_dir = local_dir

        self._model_path = f"{self._local_dir}/model.pt"
        self._adata_path = f"{self._local_dir}/adata.h5ad"
        self._adata_full_path = f"{self._local_dir}/adata_full.h5ad"

        if isinstance(model_card, ModelCard):
            self._model_card = model_card
        elif isinstance(model_card, str):
            content = Path(model_card).read_text()
            self._model_card = ModelCard(content)
        else:
            raise TypeError("`model_card` data type not understood")

        # lazy load - these are not loaded until accessed
        self._model = None
        self._adata = None
        self._adata_full = None

    def push_to_huggingface_hub(
        self, repo_name: str, repo_token_path: str, repo_create: bool
    ):
        """Placeholder docstring. TODO complete"""
        repo_token = Path(repo_token_path).read_text()
        if repo_create:
            create_repo(repo_name, token=repo_token)
        api = HfApi()
        api.upload_file(
            path_or_fileobj=self._model_path,
            path_in_repo=self._model_path.split("/")[-1],
            repo_id=repo_name,
            token=repo_token,
        )
        api.upload_file(
            path_or_fileobj=self._adata_path,
            path_in_repo=self._adata_path.split("/")[-1],
            repo_id=repo_name,
            token=repo_token,
        )
        self._model_card.push_to_hub(repo_name, token=repo_token)

    @classmethod
    def pull_from_huggingface_hub(cls, repo_id: str):
        """Placeholder docstring. TODO complete"""
        cache_dir = snapshot_download(
            repo_id=repo_id, allow_patterns=["model.pt", "adata.h5ad"]
        )
        model_card = ModelCard.load(repo_id)
        return cls(cache_dir, model_card)

    def __repr__(self):
        def eval_obj(obj):
            return "No" if obj is None else "Yes"

        return (
            "HubModel with:\n"
            f"local_dir: {self._local_dir}\n"
            f"model loaded? {eval_obj(self._model)}\n"
            f"adata loaded? {eval_obj(self._adata)}\n"
            f"adata_full loaded? {eval_obj(self._adata_full)}\n"
            f"model card:\n{self.model_card}"
        )

    @property
    def model_card(self) -> ModelCard:
        """Placeholder docstring. TODO complete"""
        return self._model_card

    @property
    def model(self) -> Type[BaseModelClass]:
        """Placeholder docstring. TODO complete"""
        if self._model is None:
            self.load_model()
        return self._model

    @property
    def adata(self) -> AnnData:
        """Placeholder docstring. TODO complete"""
        if self._adata is None:
            self.read_adata()
        return self._adata

    @property
    def adata_full(self) -> Optional[AnnData]:
        """Placeholder docstring. TODO complete"""
        if self._adata_full is None:
            self.read_adata_full()
        return self._adata_full

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
        # -> need to download it first
        raise NotImplementedError()


def list_all_models(
    do_print: bool = True, print_detailed: bool = False
) -> List[ModelInfo]:
    """Placeholder docstring. TODO complete"""
    filt = ModelFilter(library=HF_LIBRARY_NAME)
    api = HfApi()
    all_models = api.list_models(filter=filt)
    if do_print:
        for m in all_models:
            print(m if print_detailed else str(m))
    return all_models
