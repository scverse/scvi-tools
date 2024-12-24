from __future__ import annotations

import importlib
import json
import logging
import os
import tempfile
import warnings
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING

import anndata
import mudata
import rich
from huggingface_hub import ModelCard, snapshot_download
from rich.markdown import Markdown

from scvi import settings
from scvi.data import cellxgene
from scvi.data._download import _download
from scvi.hub._metadata import HubMetadata, HubModelCardHelper
from scvi.utils import dependencies

from ._constants import _SCVI_HUB

if TYPE_CHECKING:
    from anndata import AnnData

    from scvi.model.base import BaseModelClass

logger = logging.getLogger(__name__)


class HubModel:
    """Wrapper for :class:`~scvi.model.base.BaseModelClass` backed by HuggingFace Hub.

    Parameters
    ----------
    repo_name
        ID of the huggingface repo where this model is uploaded
    local_dir
        Local directory where the data and pre-trained model reside.
    metadata
        Either an instance of :class:`~scvi.hub.HubMetadata` that contains the required metadata
        for this model, or a path to a file on disk where this metadata can be read from.
    model_card
        The model card for this pre-trained model. Model card is a markdown file that describes the
        pre-trained model/data and is displayed on HuggingFace. This can be either an instance of
        :class:`~huggingface_hub.ModelCard` or an instance of :class:`~scvi.hub.HubModelCardHelper`
        that wraps the model card or a path to a file on disk where the model card can be read
        from.

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/hub/scvi_hub_intro_and_download`
    2. :doc:`/tutorials/notebooks/hub/scvi_hub_upload_and_large_files`
    """

    def __init__(
        self,
        local_dir: str,
        repo_name: str | None = None,
        metadata: HubMetadata | str | None = None,
        model_card: HubModelCardHelper | ModelCard | str | None = None,
    ):
        self._local_dir = local_dir
        self._repo_name = repo_name

        self._model_path = f"{self._local_dir}/model.pt"
        self._adata_path = f"{self._local_dir}/adata.h5ad"
        self._mudata_path = f"{self._local_dir}/mdata.h5mu"
        self._large_training_adata_path = f"{self._local_dir}/large_training_adata.h5ad"
        self._large_training_mudata_path = f"{self._local_dir}/large_training_mudata.h5mu"

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

    def save(self, overwrite: bool = False) -> None:
        """Save the model card and metadata to the model directory.

        Parameters
        ----------
        overwrite
            Whether to overwrite existing files.
        """
        card_path = os.path.join(self._local_dir, _SCVI_HUB.MODEL_CARD_FILE_NAME)
        if os.path.isfile(card_path) and not overwrite:
            raise FileExistsError(
                f"Model card already exists at {card_path}. To overwrite, pass `overwrite=True`."
            )
        self.model_card.save(card_path)

        metadata_path = os.path.join(self._local_dir, _SCVI_HUB.METADATA_FILE_NAME)
        self.metadata.save(metadata_path, overwrite=overwrite)

    @dependencies("huggingface_hub")
    def push_to_huggingface_hub(
        self,
        repo_name: str,
        repo_token: str | None = None,
        repo_create: bool = False,
        repo_create_kwargs: dict | None = None,
        collection_name: str | None = None,
        push_anndata: bool = True,
        **kwargs,
    ):
        """Push this model to huggingface.

        If the dataset is too large to upload to huggingface, this will raise an
        exception prompting the user to upload the data elsewhere. Otherwise, the
        data, model card, and metadata are all uploaded to the given model repo.

        Parameters
        ----------
        repo_name
            ID of the huggingface repo where this model needs to be uploaded
        repo_token
            huggingface API token with write permissions if None uses token in HfFolder.get_token()
        repo_create
            Whether to create the repo
        repo_create_kwargs
            Keyword arguments passed into :meth:`~huggingface_hub.create_repo` if
            ``repo_create=True``.
        collection_name
            The name of the collection to which the model belongs.
        push_anndata
            Whether to push the :class:`~anndata.AnnData` object associated with the model.
        **kwargs
            Additional keyword arguments passed into :meth:`~huggingface_hub.HfApi.upload_file`.
        """
        from huggingface_hub import HfApi, HfFolder, add_collection_item, create_repo

        self._repo_name = repo_name

        if os.path.isfile(self._adata_path) and (
            os.path.getsize(self._adata_path) >= _SCVI_HUB.MAX_HF_UPLOAD_SIZE
        ):
            raise ValueError(
                "Dataset is too large to upload to the Model. \
                Please refer to scvi-tools tutorials for how to handle this case."
            )
        if os.path.isfile(self._mudata_path) and (
            os.path.getsize(self._mudata_path) >= _SCVI_HUB.MAX_HF_UPLOAD_SIZE
        ):
            raise ValueError(
                "Dataset is too large to upload to the Model. \
                Please refer to scvi-tools tutorials for how to handle this case."
            )
        if repo_token is None:
            repo_token = HfFolder.get_token()
        elif os.path.isfile(repo_token):
            repo_token = Path(repo_token).read_text()
        if repo_create:
            repo_create_kwargs = repo_create_kwargs or {}
            create_repo(repo_name, token=repo_token, **repo_create_kwargs)
        api = HfApi(token=repo_token)
        # upload the model card
        self.model_card.push_to_hub(repo_name, token=repo_token)
        # upload the model
        # upload the metadata
        api.upload_file(
            path_or_fileobj=json.dumps(asdict(self.metadata), indent=4).encode(),
            path_in_repo=_SCVI_HUB.METADATA_FILE_NAME,
            repo_id=repo_name,
            token=repo_token,
            **kwargs,
        )
        if not push_anndata:
            kwargs["ignore_patterns"] = ["*.h5ad", "*.h5mu"]
            api.upload_folder(
                folder_path=self._local_dir,
                repo_id=repo_name,
                token=repo_token,
                **kwargs,
            )

        if collection_name == "test":
            collection_slug = "scvi-tools/test-674f56b9eb86e62d57eac5cf"
        elif "SCANVI" in self.metadata.model_cls_name:
            collection_slug = "scvi-tools/scanvi-673c3a4aabddf849496e9079"
        elif "CondSCVI" in self.metadata.model_cls_name:
            collection_slug = "scvi-tools/destvi-673c3dbf537347953810a215"
        elif "SCVI" in self.metadata.model_cls_name:
            collection_slug = "scvi-tools/scvi-673c2c0f2bf4163ef14d018d"
        elif "TOTALVI" in self.metadata.model_cls_name:
            collection_slug = "scvi-tools/totalvi-673c3d67e2882005a1d180c1"
        elif "Stereoscope" in self.metadata.model_cls_name:
            collection_slug = "scvi-tools/stereoscope-673c3ddcf1f9f7542b8819d6"
        else:
            warnings.warn(
                "No collection found for this model."
                f"Please request a new collection for {self.metadata.model_cls_name}.",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )
            collection_slug = None

        if collection_slug is not None:
            add_collection_item(
                collection_slug=collection_slug,
                item_id=repo_name,
                item_type="model",
                exists_ok=True,
                token=repo_token,
            )

    @classmethod
    def pull_from_huggingface_hub(
        cls,
        repo_name: str,
        cache_dir: str | None = None,
        revision: str | None = None,
        pull_anndata: bool = True,
        **kwargs,
    ):
        """Download the given model repo from huggingface.

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
            The revision to pull from the repo. This can be a branch name, a tag, or a full-length
            commit hash. If None, the default (latest) revision is pulled.
        pull_anndata
            Whether to pull the :class:`~anndata.AnnData` object associated with the model. If
            ``True`` but the file does not exist, will fail silently.
        kwargs
            Additional keyword arguments to pass to :meth:`~huggingface_hub.snapshot_download`.
        """
        if revision is None:
            warnings.warn(
                "No revision was passed, so the default (latest) revision will be used.",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )
        filenames = ["model.pt", _SCVI_HUB.METADATA_FILE_NAME]
        if pull_anndata:
            filenames.append("adata.h5ad")
            filenames.append("mdata.h5mu")

        snapshot_folder = snapshot_download(
            repo_id=repo_name,
            allow_patterns=filenames,
            cache_dir=cache_dir,
            revision=revision,
            **kwargs,
        )
        model_card = ModelCard.load(repo_name)
        return cls(snapshot_folder, model_card=model_card, repo_name=repo_name)

    @dependencies("boto3")
    def push_to_s3(
        self,
        s3_bucket: str,
        s3_path: str,
        push_anndata: bool = True,
        **kwargs,
    ):
        """Upload the :class:`~scvi.hub.HubModel` to an S3 bucket.

        Requires `boto3 <https://boto3.amazonaws.com/v1/documentation/api/latest/index.html>`_ to
        be installed.

        Parameters
        ----------
        s3_bucket
            The S3 bucket to which to upload the model.
        s3_path
            The S3 path where the model will be saved.
        push_anndata
            Whether to push the :class:`~anndata.AnnData` object associated with the model.
        **kwargs
            Keyword arguments passed into :func:`~boto3.client`.
        """
        from boto3 import client

        self.save(overwrite=True)
        s3 = client("s3", **kwargs)

        card_local_path = os.path.join(self._local_dir, _SCVI_HUB.MODEL_CARD_FILE_NAME)
        card_s3_path = os.path.join(s3_path, _SCVI_HUB.MODEL_CARD_FILE_NAME)
        s3.upload_file(card_local_path, s3_bucket, card_s3_path)

        metadata_local_path = os.path.join(self._local_dir, _SCVI_HUB.METADATA_FILE_NAME)
        metadata_s3_path = os.path.join(s3_path, _SCVI_HUB.METADATA_FILE_NAME)
        s3.upload_file(metadata_local_path, s3_bucket, metadata_s3_path)

        model_s3_path = os.path.join(s3_path, "model.pt")
        s3.upload_file(self._model_path, s3_bucket, model_s3_path)

        if push_anndata:
            if not os.path.isfile(self._adata_path) and not os.path.isfile(self._mudata_path):
                raise ValueError(
                    f"No AnnData file found at {self._adata_path} or {self._mudata_path}. "
                    "Please provide an AnnData file or set `push_anndata=False`."
                )
            if os.path.isfile(self._adata_path):
                adata_s3_path = os.path.join(s3_path, "adata.h5ad")
                s3.upload_file(self._adata_path, s3_bucket, adata_s3_path)
            else:
                mudata_s3_path = os.path.join(s3_path, "mudata.h5mu")
                s3.upload_file(self._mudata_path, s3_bucket, mudata_s3_path)

    @classmethod
    @dependencies("boto3")
    def pull_from_s3(
        cls,
        s3_bucket: str,
        s3_path: str,
        pull_anndata: bool = True,
        cache_dir: str | None = None,
        unsigned: bool = False,
        **kwargs,
    ) -> HubModel:
        """Download a :class:`~scvi.hub.HubModel` from an S3 bucket.

        Requires `boto3 <https://boto3.amazonaws.com/v1/documentation/api/latest/index.html>`_ to
        be installed.

        Parameters
        ----------
        s3_bucket
            The S3 bucket from which to download the model.
        s3_path
            The S3 path to the saved model.
        pull_anndata
            Whether to pull the :class:`~anndata.AnnData` object associated with the model.
        cache_dir
            The directory where the downloaded model files will be cached. Defaults to a temporary
            directory created with :func:`tempfile.mkdtemp`.
        unsigned
            Whether to use unsigned requests. If ``True`` and ``config`` is passed in ``kwargs``,
            ``config`` will be overwritten.
        **kwargs
            Keyword arguments passed into :func:`~boto3.client`.

        Returns
        -------
        The pretrained model specified by the given S3 bucket and path.
        """
        from boto3 import client
        from botocore import UNSIGNED, config

        if unsigned:
            kwargs["config"] = config.Config(signature_version=UNSIGNED)
        cache_dir = cache_dir or tempfile.mkdtemp()
        s3 = client("s3", **kwargs)

        card_s3_path = os.path.join(s3_path, _SCVI_HUB.MODEL_CARD_FILE_NAME)
        card_local_path = os.path.join(cache_dir, _SCVI_HUB.MODEL_CARD_FILE_NAME)
        s3.download_file(s3_bucket, card_s3_path, card_local_path)

        metadata_s3_path = os.path.join(s3_path, _SCVI_HUB.METADATA_FILE_NAME)
        metadata_local_path = os.path.join(cache_dir, _SCVI_HUB.METADATA_FILE_NAME)
        s3.download_file(s3_bucket, metadata_s3_path, metadata_local_path)

        model_s3_path = os.path.join(s3_path, "model.pt")
        model_local_path = os.path.join(cache_dir, "model.pt")
        s3.download_file(s3_bucket, model_s3_path, model_local_path)

        if pull_anndata:
            try:
                adata_s3_path = os.path.join(s3_path, "adata.h5ad")
                adata_local_path = os.path.join(cache_dir, "adata.h5ad")
                s3.download_file(s3_bucket, adata_s3_path, adata_local_path)
            except s3.exceptions.ClientError:
                mudata_s3_path = os.path.join(s3_path, "mudata.h5mu")
                mudata_local_path = os.path.join(cache_dir, "mdata.h5mu")
                s3.download_file(s3_bucket, mudata_s3_path, mudata_local_path)

        model_card = ModelCard.load(card_local_path)
        return cls(cache_dir, model_card=model_card)

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
    def local_dir(self) -> str:
        """The local directory where the data and pre-trained model reside."""
        return self._local_dir

    @property
    def repo_name(self) -> str:
        """The local directory where the data and pre-trained model reside."""
        return self._repo_name

    @property
    def metadata(self) -> HubMetadata:
        """The metadata for this model."""
        return self._metadata

    @property
    def model_card(self) -> ModelCard:
        """The model card for this model."""
        return self._model_card

    @property
    def model(self) -> type[BaseModelClass]:
        """Returns the model object for this hub model.

        If the model has not been loaded yet, this will call :meth:`~scvi.hub.HubModel.load_model`.
        Otherwise, it will simply return the loaded model.
        """
        if self._model is None:
            self.load_model()
        return self._model

    @property
    def adata(self) -> AnnData | None:
        """Returns the data for this model.

        If the data has not been loaded yet, this will call :meth:`~scvi.hub.HubModel.read_adata`.
        Otherwise, it will simply return the loaded data.
        """
        if self._adata is None:
            self.read_adata()
        return self._adata

    @property
    def large_training_adata(self) -> AnnData | None:
        """Returns the training data for this model.

        If the data has not been loaded yet, this will call
        :meth:`~scvi.hub.HubModel.read_large_training_adata`, which will attempt to download from
        the source url. Otherwise, it will simply return the loaded data.
        """
        if self._large_training_adata is None:
            self.read_large_training_adata()
        return self._large_training_adata

    def load_model(
        self,
        adata: AnnData | None = None,
        accelerator: str = "auto",
        device: str | int | None = "auto",
    ):
        """Loads the model.

        Parameters
        ----------
        adata
            The data to load the model with, if not None. If None, we'll try to load the model
            using the data at ``self._adata_path``. If that file does not exist, we'll try to load
            the model using :meth:`~scvi.hub.HubModel.large_training_adata`. If that does not exist
            either, we'll error out.
        %(param_accelerator)s
        %(param_device)s
        """
        logger.info("Loading model...")
        # get the class name for this model (e.g. TOTALVI)
        model_cls_name = self.metadata.model_cls_name
        python_module = importlib.import_module(self.metadata.model_parent_module)
        model_cls = getattr(python_module, model_cls_name)
        if (
            adata is not None
            or os.path.isfile(self._adata_path)
            or os.path.isfile(self._mudata_path)
        ):
            self._model = model_cls.load(
                os.path.dirname(self._model_path),
                adata=adata,
                accelerator=accelerator,
                device=device,
            )
        else:
            # in this case, we must download the large training adata if it exists in the model
            # card; otherwise, we error out. Note that the call below faults in
            # self.large_training_adata if it is None
            if self.large_training_adata is None:
                raise ValueError(
                    "Could not find any dataset to load the model with. Either provide "
                    "a dataset on disk or a url to download the data in the model "
                    "card, or pass an `adata` to this method. See scvi-tools tutorials "
                    "for more details."
                )
            else:
                self._model = model_cls.load(
                    os.path.dirname(self._model_path),
                    adata=self.large_training_adata,
                    accelerator=accelerator,
                    device=device,
                )

    def read_adata(self) -> None:
        """Reads the data from disk (``self._adata_path``).

        Reads if it exists. Otherwise, this is a no-op.
        """
        if os.path.isfile(self._adata_path):
            logger.info("Reading adata...")
            self._adata = anndata.read_h5ad(self._adata_path)
        if os.path.isfile(self._mudata_path):
            logger.info("Reading adata...")
            self._adata = mudata.read_h5mu(self._mudata_path)
        else:
            logger.info("No data found on disk. Skipping...")

    def read_mudata(self) -> None:
        """Reads the data from disk (``self._mudata_path``).

        Reads if it exists. Otherwise, this is a no-op.
        """
        if os.path.isfile(self._mudata_path):
            logger.info("Reading adata...")
            self._adata = anndata.read_h5ad(self._mudata_path)
        else:
            logger.info("No data found on disk. Skipping...")

    def read_large_training_adata(self) -> None:
        """Downloads the large training adata.

        If it exists, then load it into memory. Otherwise, this is a no-op.

        Notes
        -----
        The large training data url can be a cellxgene explorer session url. However it cannot be a
        self-hosted session. In other words, it must be from the cellxgene portal
        (https://cellxgene.cziscience.com/).
        """
        training_data_url = self.metadata.training_data_url
        if self.metadata.model_cls_name in ["TOTALVI", "MULTIVI"]:
            large_training_adata_path = self._large_training_mudata_path
        else:
            large_training_adata_path = self._large_training_adata_path
        if training_data_url is not None:
            logger.info(
                f"Downloading large training dataset from this url:\n{training_data_url}..."
            )
            # Add multi-modal models here.
            dn = Path(large_training_adata_path).parent.as_posix()
            fn = Path(large_training_adata_path).name
            url_parts = training_data_url.split("/")
            url_last_part = url_parts[-2] if url_parts[-1] == "" else url_parts[-1]
            if url_last_part.endswith(".cxg"):
                _ = cellxgene(training_data_url, fn, dn, return_path=True)
            else:
                _download(training_data_url, dn, fn)
            logger.info("Reading large training data...")
        else:
            # Access file from DVC storage.
            from dvc.api import DVCFileSystem

            if self.metadata.model_cls_name in ["TOTALVI", "MULTIVI"]:
                suffix = "h5mu"
            else:
                suffix = "h5ad"

            fs = DVCFileSystem(
                "https://github.com/YosefLab/scvi-hub-models", rev="main", remote="s3_remote"
            )
            if self.repo_name is not None:
                fs.get_file(
                    f"data/{self.repo_name.split('/')[1].rsplit('_', 1)[0]}.{suffix}",
                    large_training_adata_path,
                )
        if os.path.exists(large_training_adata_path):
            if large_training_adata_path.endswith(".h5mu"):
                self._large_training_adata = mudata.read_h5mu(large_training_adata_path)
            else:
                self._large_training_adata = anndata.read_h5ad(large_training_adata_path)
        else:
            self._large_training_adata = None
