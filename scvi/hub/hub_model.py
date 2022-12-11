import importlib
import logging
import os
from pathlib import Path
from typing import List, Optional, Type, Union

import anndata
import rich
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
from rich.markdown import Markdown

from scvi.data._download import _download
from scvi.hub.hub_metadata import HubMetadata, HubModelCardHelper
from scvi.model.base import BaseModelClass

from ._constants import (
    HF_LIBRARY_NAME,
    MAX_HF_UPLOAD_SIZE,
    METADATA_FILE_NAME,
    MODEL_CARD_FILE_NAME,
)

logger = logging.getLogger(__name__)


class HubModel:
    """Placeholder docstring. TODO complete."""

    def __init__(
        self,
        local_dir: str,
        metadata: Optional[Union[HubMetadata, str]] = None,
        model_card: Optional[Union[HubModelCardHelper, ModelCard, str]] = None,
    ):
        self._local_dir = local_dir

        self._model_path = f"{self._local_dir}/model.pt"
        self._adata_path = f"{self._local_dir}/adata.h5ad"
        self._adata_large_path = f"{self._local_dir}/adata_large.h5ad"
        self._metadata_path = f"{self._local_dir}/{METADATA_FILE_NAME}"
        self._model_card_path = f"{self._local_dir}/{MODEL_CARD_FILE_NAME}"

        # lazy load - these are not loaded until accessed
        self._model = None
        self._adata = None
        self._adata_large = None

        # get the metadata from the parameters or from the disk
        if isinstance(metadata, HubMetadata):
            self._metadata = metadata
        elif isinstance(metadata, str) or os.path.isfile(self._metadata_path):
            path = metadata if isinstance(metadata, str) else self._metadata_path
            content = Path(path).read_text()
            self._metadata = HubMetadata.from_json(content)
        else:
            raise ValueError("No metadata found")

        # get the model card from the parameters or from the disk
        if isinstance(model_card, HubModelCardHelper):
            self._model_card = model_card.model_card
        elif isinstance(model_card, ModelCard):
            self._model_card = model_card
        elif isinstance(model_card, str) or os.path.isfile(self._model_card_path):
            path = model_card if isinstance(model_card, str) else self._model_card_path
            content = Path(path).read_text()
            self._model_card = ModelCard(content)
        else:
            raise ValueError("No model card found")

    def push_to_huggingface_hub(
        self, repo_name: str, repo_token_path: str, repo_create: bool
    ):
        """Placeholder docstring. TODO complete."""
        if os.path.isfile(self._adata_path) and (
            os.path.getsize(self._adata_path) >= MAX_HF_UPLOAD_SIZE
        ):
            raise ValueError(
                "Dataset is too large to upload to the Model. \
                Please refer to scvi-tools tutorials for how to handle this case."
            )
        repo_token = Path(repo_token_path).read_text()
        if repo_create:
            create_repo(repo_name, token=repo_token)
        api = HfApi()
        # upload the model
        api.upload_file(
            path_or_fileobj=self._model_path,
            path_in_repo=self._model_path.split("/")[-1],
            repo_id=repo_name,
            token=repo_token,
        )
        # upload the data if it exists
        if os.path.isfile(self._adata_path):
            api.upload_file(
                path_or_fileobj=self._adata_path,
                path_in_repo=self._adata_path.split("/")[-1],
                repo_id=repo_name,
                token=repo_token,
            )
        # upload the metadata and model card
        api.upload_file(
            # path_or_fileobj=json.dumps(asdict(self.metadata), indent=4).encode(),
            path_or_fileobj=self.metadata.to_json(indent=4).encode(),
            path_in_repo=METADATA_FILE_NAME,
            repo_id=repo_name,
            token=repo_token,
        )
        # upload the model card
        self.model_card.push_to_hub(repo_name, token=repo_token)

    @classmethod
    def pull_from_huggingface_hub(cls, repo_name: str):
        """Placeholder docstring. TODO complete."""
        cache_dir = snapshot_download(
            repo_id=repo_name,
            allow_patterns=["model.pt", "adata.h5ad", METADATA_FILE_NAME],
        )
        model_card = ModelCard.load(repo_name)
        return cls(cache_dir, model_card)

    def __repr__(self):
        def eval_obj(obj):
            return "No" if obj is None else "Yes"

        print(
            "HubModel with:\n"
            f"local_dir: {self._local_dir}\n"
            f"model loaded? {eval_obj(self._model)}\n"
            f"adata loaded? {eval_obj(self._adata)}\n"
            f"adata_large loaded? {eval_obj(self._adata_large)}\n"
            f"metadata:\n{self.metadata}\n"
            f"model_card:"
        )
        rich.print(Markdown(self.model_card.content))
        return ""

    @property
    def metadata(self) -> HubMetadata:
        """Placeholder docstring. TODO complete."""
        return self._metadata

    @property
    def model_card(self) -> ModelCard:
        """Placeholder docstring. TODO complete."""
        return self._model_card

    @property
    def model(self) -> Type[BaseModelClass]:
        """Placeholder docstring. TODO complete."""
        if self._model is None:
            self.load_model()
        return self._model

    @property
    def adata(self) -> AnnData:
        """Placeholder docstring. TODO complete."""
        if self._adata is None:
            self.read_adata()
        return self._adata

    @property
    def adata_large(self) -> Optional[AnnData]:
        """Placeholder docstring. TODO complete."""
        if self._adata_large is None:
            self.read_adata_large()
        return self._adata_large

    def load_model(
        self,
        adata: Optional[AnnData] = None,
    ):
        """Placeholder docstring. TODO complete."""
        logger.info("Loading model...")
        # get the class name for this model (e.g. TOTALVI)
        torch_model = torch.load(self._model_path, map_location="cpu")
        cls_name = torch_model["attr_dict"]["registry_"]["model_name"]
        python_module = importlib.import_module(self.metadata.model_parent_module)
        model_cls = getattr(python_module, cls_name)
        if adata is not None or os.path.isfile(self._adata_path):
            self._model = model_cls.load(os.path.dirname(self._model_path), adata=adata)
        else:
            # in this case, we must download the large adata if it exists in the model card; otherwise, we error out
            # the call below faults in self.adata_large if it is None
            if self.adata_large is None:
                raise ValueError(
                    "Could not find any dataset to load the model with.\
                    Either provide a dataset on disk or a url to download the data in the model card.\
                    See scvi-tools tutorials for more details."
                )
            else:
                self._model = model_cls.load(
                    os.path.dirname(self._model_path),
                    adata=self.adata_large,
                )

    def read_adata(self):
        """Placeholder docstring. TODO complete."""
        if os.path.isfile(self._adata_path):
            logger.info("Reading adata...")
            self._adata = anndata.read_h5ad(self._adata_path)
        else:
            logger.info("No data found on disk. Skipping...")

    def read_adata_large(self):
        """Download the full adata, if it exists, then read it into memory."""
        large_data_url = self.metadata.large_data_url
        if large_data_url is not None:
            logger.info(
                f"Downloading large dataset from this url:\n{large_data_url}..."
            )
            dn = Path(self._adata_large_path).parent.as_posix()
            fn = Path(self._adata_large_path).name
            _download(large_data_url, dn, fn)
            logger.info("Reading large data...")
            self._adata_large = anndata.read_h5ad(self._adata_large_path)
        else:
            logger.info("No large_data_url found in the model card. Skipping...")


def list_all_models(
    do_print: bool = True, print_detailed: bool = False
) -> Optional[List[ModelInfo]]:
    """Placeholder docstring. TODO complete."""
    filt = ModelFilter(library=HF_LIBRARY_NAME)
    api = HfApi()
    all_models = api.list_models(filter=filt)
    if do_print:
        for m in all_models:
            print("\n")
            logger.info(m.__repr__() if print_detailed else m)
    return None if do_print else all_models
