import importlib
import json
import logging
import os
import warnings
from dataclasses import asdict
from pathlib import Path
from typing import Optional, Type, Union

import anndata
import rich
from anndata import AnnData
from huggingface_hub import HfApi, ModelCard, create_repo, snapshot_download
from rich.markdown import Markdown

from scvi.data._download import _download
from scvi.hub.hub_metadata import HubMetadata, HubModelCardHelper
from scvi.model.base import BaseModelClass
from scvi.model.base._utils import _load_saved_files

from ._constants import _SCVI_HUB

logger = logging.getLogger(__name__)


class HubModel:
    """
    Provides functionality to interact with the scvi-hub backed by `huggingface <https://huggingface.co/models>`_.

    Parameters
    ----------
    local_dir
        Local directory where the data and pre-trained model reside.
    metadata
        Either an instance of :class:`~scvi.hub.HubMetadata` that contains the required metadata for this model,
        or a path to a file on disk where this metadata can be read from.
    model_card
        The model card for this pre-trained model. Model card is a markdown file that describes the pre-trained
        model/data and is displayed on huggingface. This can be either an instance of
        :class:`~huggingface_hub.ModelCard` or an instance of :class:`~scvi.hub.HubModelCardHelper` that wraps
        the model card or a path to a file on disk where the model card can be read from.

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/scvi_hub_intro_and_download`
    2. :doc:`/tutorials/notebooks/scvi_hub_upload_and_large_files`
    """

    def __init__(
        self,
        local_dir: str,
        metadata: Optional[Union[HubMetadata, str]] = None,
        model_card: Optional[Union[HubModelCardHelper, ModelCard, str]] = None,
    ):
        self._local_dir = local_dir

        self._model_path = f"{self._local_dir}/model.pt"
        self._adata_path = f"{self._local_dir}/adata.h5ad"
        self._large_training_adata_path = f"{self._local_dir}/large_training_adata.h5ad"

        # lazy load - these are not loaded until accessed
        self._model = None
        self._adata = None
        self._large_training_adata = None

        # get the metadata from the parameters or from the disk
        metadata_path = f"{self._local_dir}/{_SCVI_HUB.METADATA_FILE_NAME}"
        if isinstance(metadata, HubMetadata):
            self._metadata = metadata
        elif isinstance(metadata, str) or os.path.isfile(metadata_path):
            path = metadata if isinstance(metadata, str) else metadata_path
            content = Path(path).read_text()
            content_dict = json.loads(content)
            self._metadata = HubMetadata(**content_dict)
        else:
            raise ValueError("No metadata found")

        # get the model card from the parameters or from the disk
        model_card_path = f"{self._local_dir}/{_SCVI_HUB.MODEL_CARD_FILE_NAME}"
        if isinstance(model_card, HubModelCardHelper):
            self._model_card = model_card.model_card
        elif isinstance(model_card, ModelCard):
            self._model_card = model_card
        elif isinstance(model_card, str) or os.path.isfile(model_card_path):
            path = model_card if isinstance(model_card, str) else model_card_path
            content = Path(path).read_text()
            self._model_card = ModelCard(content)
        else:
            raise ValueError("No model card found")

    def push_to_huggingface_hub(
        self, repo_name: str, repo_token: str, repo_create: bool
    ):
        """
        Push this model to huggingface.

        If the dataset is too large to upload to huggingface, this will raise an
        exception prompting the user to upload the data elsewhere. Otherwise, the
        data, model card, and metadata are all uploaded to the given model repo.

        Parameters
        ----------
        repo_name
            ID of the huggingface repo where this model needs to be uploaded
        repo_token
            huggingface API token with write permissions
        repo_create
            Whether to create the repo
        """
        if os.path.isfile(self._adata_path) and (
            os.path.getsize(self._adata_path) >= _SCVI_HUB.MAX_HF_UPLOAD_SIZE
        ):
            raise ValueError(
                "Dataset is too large to upload to the Model. \
                Please refer to scvi-tools tutorials for how to handle this case."
            )
        if os.path.isfile(repo_token):
            repo_token = Path(repo_token).read_text()
        if repo_create:
            create_repo(repo_name, token=repo_token)
        api = HfApi()
        # upload the model card
        self.model_card.push_to_hub(repo_name, token=repo_token)
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
        # upload the metadata
        api.upload_file(
            path_or_fileobj=json.dumps(asdict(self.metadata), indent=4).encode(),
            path_in_repo=_SCVI_HUB.METADATA_FILE_NAME,
            repo_id=repo_name,
            token=repo_token,
        )

    @classmethod
    def pull_from_huggingface_hub(
        cls,
        repo_name: str,
        cache_dir: Optional[str] = None,
        revision: Optional[str] = None,
        **kwargs,
    ):
        """
        Download the given model repo from huggingface.

        The model, its card, data, metadata are downloaded to a cached location on disk
        selected by huggingface and an instance of this class is created with that info
        and returned.

        Parameters
        ----------
        repo_name
            ID of the huggingface repo where this model needs to be uploaded
        cache_dir
            The directory where the downloaded model artifacts will be cached
        revision
            The revision to pull from the repo. This can be a branch name, a tag, or a full-length commit hash.
            If None, the default (latest) revision is pulled.
        kwargs
            Additional keyword arguments to pass to :meth:`~huggingface_hub.snapshot_download`.
        """
        if revision is None:
            warnings.warn(
                "No revision was passed, so the default (latest) revision will be used."
            )
        snapshot_folder = snapshot_download(
            repo_id=repo_name,
            allow_patterns=["model.pt", "adata.h5ad", _SCVI_HUB.METADATA_FILE_NAME],
            cache_dir=cache_dir,
            revision=revision,
            **kwargs,
        )
        model_card = ModelCard.load(repo_name)
        return cls(snapshot_folder, model_card=model_card)

    def __repr__(self):
        def eval_obj(obj):
            return "No" if obj is None else "Yes"

        print(
            "HubModel with:\n"
            f"local_dir: {self._local_dir}\n"
            f"model loaded? {eval_obj(self._model)}\n"
            f"adata loaded? {eval_obj(self._adata)}\n"
            f"large_training_adata loaded? {eval_obj(self._large_training_adata)}\n"
            f"metadata:\n{self.metadata}\n"
            f"model_card:"
        )
        # TODO figure out how to print tables in rich Markdown without the replace trick
        rich.print(Markdown(self.model_card.content.replace("\n", "\n\n")))
        return ""

    @property
    def metadata(self) -> HubMetadata:
        """The metadata for this model."""
        return self._metadata

    @property
    def model_card(self) -> ModelCard:
        """The model card for this model."""
        return self._model_card

    @property
    def model(self) -> Type[BaseModelClass]:
        """
        Returns the model object for this hub model.

        If the model has not been loaded yet, this will call :meth:`~scvi.hub.HubModel.load_model`.
        Otherwise, it will simply return the loaded model.
        """
        if self._model is None:
            self.load_model()
        return self._model

    @property
    def adata(self) -> Optional[AnnData]:
        """
        Returns the data for this model.

        If the data has not been loaded yet, this will call :meth:`~scvi.hub.HubModel.read_adata`.
        Otherwise, it will simply return the loaded data.
        """
        if self._adata is None:
            self.read_adata()
        return self._adata

    @property
    def large_training_adata(self) -> Optional[AnnData]:
        """
        Returns the training data for this model, which might be too large to reside within the hub model.

        If the data has not been loaded yet, this will call :meth:`~scvi.hub.HubModel.read_large_training_adata`,
        which will attempt to download from the source url. Otherwise, it will simply return the loaded data.
        """
        if self._large_training_adata is None:
            self.read_large_training_adata()
        return self._large_training_adata

    def load_model(
        self,
        adata: Optional[AnnData] = None,
    ):
        """
        Loads the model.

        Parameters
        ----------
        adata
            The data to  load the model with, if not None. If None, we'll try to load the model using the data
            at ``self._adata_path``. If that file does not exist, we'll try to load the model using
            :meth:`~scvi.hub.HubModel.large_training_adata`. If that does not exist either, we'll error out.
        """
        logger.info("Loading model...")
        # get the class name for this model (e.g. TOTALVI)
        attr_dict, _, _, _ = _load_saved_files(self._local_dir, load_adata=False)
        cls_name = attr_dict["registry_"]["model_name"]
        python_module = importlib.import_module(self.metadata.model_parent_module)
        model_cls = getattr(python_module, cls_name)
        if adata is not None or os.path.isfile(self._adata_path):
            self._model = model_cls.load(os.path.dirname(self._model_path), adata=adata)
        else:
            # in this case, we must download the large training adata if it exists in the model card; otherwise,
            # we error out. Note that the call below faults in self.large_training_adata if it is None
            if self.large_training_adata is None:
                raise ValueError(
                    "Could not find any dataset to load the model with.\
                    Either provide a dataset on disk or a url to download the data in the model card,\
                    or pass an adata to this method.\
                    See scvi-tools tutorials for more details."
                )
            else:
                self._model = model_cls.load(
                    os.path.dirname(self._model_path),
                    adata=self.large_training_adata,
                )

    def read_adata(self):
        """Reads the data from disk (``self._adata_path``) if it exists. Otherwise, this is a no-op."""
        if os.path.isfile(self._adata_path):
            logger.info("Reading adata...")
            self._adata = anndata.read_h5ad(self._adata_path)
        else:
            logger.info("No data found on disk. Skipping...")

    def read_large_training_adata(self):
        """Downloads the large training adata, if it exists, then load it into memory. Otherwise, this is a no-op."""
        training_data_url = self.metadata.training_data_url
        if training_data_url is not None:
            logger.info(
                f"Downloading large training dataset from this url:\n{training_data_url}..."
            )
            dn = Path(self._large_training_adata_path).parent.as_posix()
            fn = Path(self._large_training_adata_path).name
            _download(training_data_url, dn, fn)
            logger.info("Reading large training data...")
            self._large_training_adata = anndata.read_h5ad(
                self._large_training_adata_path
            )
        else:
            logger.info("No training_data_url found in the model card. Skipping...")
