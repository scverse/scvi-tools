import json
from dataclasses import dataclass
from typing import List, Optional

import torch
from dataclasses_json import dataclass_json
from huggingface_hub import ModelCard, ModelCardData

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
from ._utils import _get_counts_and_latent_info
from .model_card_template import template


@dataclass_json
@dataclass
class HubMetadata:
    """Placeholder docstring. TODO complete."""

    data_cell_count: int
    data_gene_count: int
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
        data_cell_count: Optional[int] = None,
        data_gene_count: Optional[int] = None,
        **kwargs,
    ):
        """Placeholder docstring. TODO complete."""
        cell_count, gene_count = _get_counts_and_latent_info(
            local_dir, data_cell_count, data_gene_count
        )
        torch_model = torch.load(f"{local_dir}/model.pt")
        scvi_version = torch_model["attr_dict"]["registry_"]["scvi_version"]

        return cls(
            cell_count,
            gene_count,
            scvi_version,
            anndata_version,
            **kwargs,
        )


@dataclass
class HubModelCardHelper:
    """Placeholder docstring. TODO complete."""

    license_info: str
    data_cell_count: int
    data_gene_count: int
    model_cls_name: str
    model_init_params: dict
    model_setup_anndata_args: dict
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
        data_cell_count: Optional[int] = None,
        data_gene_count: Optional[int] = None,
        data_is_latent: Optional[bool] = None,
        **kwargs,
    ):
        """Placeholder docstring. TODO complete."""
        cell_count, gene_count, is_latent = _get_counts_and_latent_info(
            local_dir, data_cell_count, data_gene_count, data_is_latent
        )

        torch_model = torch.load(f"{local_dir}/model.pt")
        attr_dict = torch_model["attr_dict"]
        model_init_params = attr_dict["init_params_"]
        model_cls_name = attr_dict["registry_"]["model_name"]
        model_setup_anndata_args = attr_dict["registry_"]["setup_args"]
        scvi_version = attr_dict["registry_"]["scvi_version"]

        return cls(
            license_info,
            cell_count,
            gene_count,
            model_cls_name,
            model_init_params,
            model_setup_anndata_args,
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
            cell_count=self.data_cell_count,
            gene_count=self.data_gene_count,
            model_init_params=json.dumps(self.model_init_params, indent=4),
            model_setup_anndata_args=json.dumps(
                self.model_setup_anndata_args, indent=4
            ),
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
