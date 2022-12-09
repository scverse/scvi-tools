import json
import os
from dataclasses import dataclass
from typing import List, Optional

import anndata
import torch
from dataclasses_json import dataclass_json
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


@dataclass_json
@dataclass
class HubMetadata:
    """Placeholder docstring. TODO complete."""

    scvi_version: str
    anndata_version: str
    large_data_url: Optional[str] = None
    model_parent_module: str = DEFAULT_PARENT_MODULE

    # TODO consider checking large_data_url is valid url

    @classmethod
    def from_dir(
        cls,
        local_dir: str,
        anndata_version: str,
        **kwargs,
    ):
        """Placeholder docstring. TODO complete."""
        torch_model = torch.load(f"{local_dir}/model.pt")
        scvi_version = torch_model["attr_dict"]["registry_"]["scvi_version"]

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
    model_anndata_setup_view: str
    scvi_version: str
    anndata_version: str
    data_modalities: Optional[List[str]] = None
    tissues: Optional[List[str]] = None
    data_is_annotated: Optional[bool] = None
    data_is_latent: Optional[bool] = None
    large_data_url: Optional[str] = None
    model_parent_module: str = DEFAULT_PARENT_MODULE
    description: str = DEFAULT_MISSING_FIELD
    references: str = DEFAULT_MISSING_FIELD

    def __post_init__(self):
        self.data_modalities = self.data_modalities or []
        self.tissues = self.tissues or []
        self.model_card = self._to_model_card()

    @classmethod
    def from_dir(
        cls,
        local_dir: str,
        license_info: str,
        anndata_version: str,
        torch_map_location: str = "cpu",
        data_is_latent: Optional[bool] = None,
        **kwargs,
    ):
        """Placeholder docstring. TODO complete."""
        torch_model = torch.load(
            f"{local_dir}/model.pt", map_location=torch_map_location
        )
        model_init_params = torch_model["attr_dict"]["init_params_"]
        registry = torch_model["attr_dict"]["registry_"]
        model_anndata_setup_view = AnnDataManager.view_registry_from_dict(
            registry, as_str=True
        )
        model_cls_name = registry["model_name"]
        scvi_version = registry["scvi_version"]

        # get `is_latent` from the param if it is given, else from adata if it on disk, else set it to None
        is_latent = data_is_latent
        if is_latent is None and os.path.isfile(f"{local_dir}/adata.h5ad"):
            adata = anndata.read_h5ad(f"{local_dir}/adata.h5ad", backed=True)
            is_latent = _is_latent(adata)

        return cls(
            license_info,
            model_cls_name,
            model_init_params,
            model_anndata_setup_view,
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
            model_anndata_setup_view=self.model_anndata_setup_view,
            model_parent_module=self.model_parent_module,
            data_is_latent=DEFAULT_MISSING_FIELD
            if self.data_is_latent is None
            else self.data_is_latent,
            large_data_url=self.large_data_url or DEFAULT_NA_FIELD,
            references=self.references,
        )

        # finally create and return the actual card
        # TODO run card.validate()? is it slow?
        return ModelCard(content)
