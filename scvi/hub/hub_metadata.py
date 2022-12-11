import json
import os
from dataclasses import dataclass, field
from typing import List, Optional

import anndata
import torch
from huggingface_hub import ModelCard, ModelCardData

from scvi.data import AnnDataManager
from scvi.data._utils import _is_latent

from ._constants import (
    ANNDATA_VERSION_TAG,
    ANNOTATED_TAG,
    DEFAULT_MISSING_FIELD,
    DEFAULT_NA_FIELD,
    DEFAULT_PARENT_MODULE,
    HF_LIBRARY_NAME,
    MODALITY_TAG,
    MODEL_CLS_NAME_TAG,
    SCVI_VERSION_TAG,
    TISSUE_TAG,
)
from .model_card_template import template


@dataclass
class HubMetadata:
    """Placeholder docstring. TODO complete."""

    scvi_version: str
    anndata_version: str
    large_data_url: Optional[str] = None
    model_parent_module: str = DEFAULT_PARENT_MODULE

    @classmethod
    def from_dict(cls, d: dict):
        """Placeholder docstring. TODO complete."""
        return cls(**d)

    @classmethod
    def from_dir(
        cls,
        local_dir: str,
        anndata_version: str,
        **kwargs,
    ):
        """Placeholder docstring. TODO complete."""
        model = torch.load(f"{local_dir}/model.pt", map_location="cpu")
        scvi_version = model["attr_dict"]["registry_"]["scvi_version"]

        return cls(
            scvi_version,
            anndata_version,
            **kwargs,
        )


@dataclass
class HubModelCardHelper:
    """Placeholder docstring. TODO complete."""

    license_info: str
    model_cls_name: str
    model_init_params: dict
    model_setup_anndata_args: dict
    model_summary_stats: dict
    model_data_registry: dict
    scvi_version: str
    anndata_version: str
    data_modalities: List[str] = field(default_factory=list)
    tissues: List[str] = field(default_factory=list)
    data_is_annotated: Optional[bool] = None
    data_is_latent: Optional[bool] = None
    large_data_url: Optional[str] = None
    model_parent_module: str = DEFAULT_PARENT_MODULE
    description: str = DEFAULT_MISSING_FIELD
    references: str = DEFAULT_MISSING_FIELD

    def __post_init__(self):
        self.model_card = self._to_model_card()

    @classmethod
    def from_dir(
        cls,
        local_dir: str,
        license_info: str,
        anndata_version: str,
        data_is_latent: Optional[bool] = None,
        **kwargs,
    ):
        """Placeholder docstring. TODO complete."""
        model = torch.load(f"{local_dir}/model.pt", map_location="cpu")
        model_init_params = model["attr_dict"]["init_params_"]
        registry = model["attr_dict"]["registry_"]
        model_cls_name = registry["model_name"]
        scvi_version = registry["scvi_version"]
        model_setup_anndata_args = registry["setup_args"]
        model_summary_stats = dict(
            AnnDataManager._get_summary_stats_from_registry(registry)
        )
        model_data_registry = dict(
            AnnDataManager._get_data_registry_from_registry(registry)
        )

        # get `is_latent` from the param if it is given, else from adata if it on disk, else set it to None
        is_latent = data_is_latent
        if is_latent is None and os.path.isfile(f"{local_dir}/adata.h5ad"):
            adata = anndata.read_h5ad(f"{local_dir}/adata.h5ad", backed=True)
            is_latent = _is_latent(adata)

        return cls(
            license_info,
            model_cls_name,
            model_init_params,
            model_setup_anndata_args,
            model_summary_stats,
            model_data_registry,
            scvi_version,
            anndata_version,
            data_is_latent=is_latent,
            **kwargs,
        )

    def _to_model_card(self) -> ModelCard:
        """Placeholder docstring. TODO complete."""
        # define tags
        tags = [
            MODEL_CLS_NAME_TAG.format(self.model_cls_name),
            SCVI_VERSION_TAG.format(self.scvi_version),
            ANNDATA_VERSION_TAG.format(self.anndata_version),
        ]
        for m in self.data_modalities:
            tags.append(MODALITY_TAG.format(m))
        for t in self.tissues:
            tags.append(TISSUE_TAG.format(t))
        if self.data_is_annotated is not None:
            tags.append(ANNOTATED_TAG.format(self.data_is_annotated))

        # define the card data, which is the header
        card_data = ModelCardData(
            license=self.license_info,
            library_name=HF_LIBRARY_NAME,
            tags=tags,
        )

        # create the content from the template
        content = template.format(
            card_data=card_data.to_yaml(),
            description=self.description,
            model_init_params=json.dumps(self.model_init_params, indent=4),
            model_setup_anndata_args=json.dumps(
                self.model_setup_anndata_args, indent=4
            ),
            model_summary_stats=AnnDataManager._view_summary_stats(
                self.model_summary_stats, as_markdown=True
            ),
            model_data_registry=AnnDataManager._view_data_registry(
                self.model_data_registry, as_markdown=True
            ),
            model_parent_module=self.model_parent_module,
            data_is_latent=DEFAULT_MISSING_FIELD
            if self.data_is_latent is None
            else self.data_is_latent,
            large_data_url=self.large_data_url or DEFAULT_NA_FIELD,
            references=self.references,
        )

        # finally create and return the actual card
        return ModelCard(content)
