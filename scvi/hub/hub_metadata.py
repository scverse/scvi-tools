import json
import logging
from dataclasses import dataclass
from typing import List, Optional

import torch
from dataclasses_json import dataclass_json
from huggingface_hub import ModelCard, ModelCardData

from ._utils import _get_cell_gene_counts
from .model_card_template import template

logger = logging.getLogger(__name__)

# TODO move to constants.py
HF_LIBRARY_NAME = "scvi-tools"
# defaults
DEFAULT_MISSING_FIELD = "To be added..."
DEFAULT_NA_FIELD = "N/A"
DEFAULT_PARENT_MODULE = "scvi.model"
# model card tags
MODEL_CLS_NAME_TAG = "model_cls_name:{}"
SCVI_VERSION_TAG = "scvi_version:{}"
ANNDATA_VERSION_TAG = "anndata_version:{}"
MODALITY_TAG = "modality:{}"
TISSUE_TAG = "tissue:{}"
ANNOTATED_TAG = "annotated:{}"


# TODO move to own file
@dataclass_json
@dataclass
class HubMetadata:
    """
    Placeholder docstring. TODO complete

    Parameters
    ----------
    data_cell_count
        number of cells in the dataset
    data_gene_count
        number of genes in the dataset
    """

    data_cell_count: int
    data_gene_count: int
    scvi_version: str
    anndata_version: str
    large_data_url: Optional[str] = None
    model_parent_module: str = DEFAULT_PARENT_MODULE

    @classmethod
    def from_dir(
        cls,
        local_dir: str,
        anndata_version: str,
        data_cell_count: Optional[int] = None,
        data_gene_count: Optional[int] = None,
        **kwargs,
    ):
        """Placeholder docstring. TODO complete"""
        cell_count, gene_count = _get_cell_gene_counts(
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


class HubModelCardHelper:
    """Placeholder docstring. TODO complete"""

    def __init__(
        self,
        license_info: str,
        data_cell_count: int,
        data_gene_count: int,
        model_cls_name: str,
        model_init_params: dict,
        model_setup_anndata_args: dict,
        scvi_version: str,
        anndata_version: str,
        data_modalities: Optional[List[str]] = None,
        tissues: Optional[List[str]] = None,
        data_is_annotated: Optional[bool] = None,
        large_data_url: Optional[str] = None,
        model_parent_module: str = DEFAULT_PARENT_MODULE,
        description: str = DEFAULT_MISSING_FIELD,
        references: str = DEFAULT_MISSING_FIELD,
    ):
        self._data_cell_count = data_cell_count
        self._data_gene_count = data_gene_count
        self._data_modalities = data_modalities or []
        self._tissues = tissues or []
        self._data_is_annotated = data_is_annotated
        self._large_data_url = large_data_url

        self._model_cls_name = model_cls_name
        self._model_init_params = model_init_params
        self._model_setup_anndata_args = model_setup_anndata_args
        self._model_parent_module = model_parent_module

        self._license_info = license_info
        self._scvi_version = scvi_version
        self._anndata_version = anndata_version
        self._description = description
        self._references = references

        self._model_card = self._to_model_card()

    @classmethod
    def from_dir(
        cls,
        local_dir: str,
        license_info: str,
        anndata_version: str,
        data_cell_count: Optional[int] = None,
        data_gene_count: Optional[int] = None,
        **kwargs,
    ):
        """Placeholder docstring. TODO complete"""
        cell_count, gene_count = _get_cell_gene_counts(
            local_dir, data_cell_count, data_gene_count
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
            **kwargs,
        )

    def _to_model_card(self) -> ModelCard:
        """Placeholder docstring. TODO complete"""
        # define tags
        tags = [
            MODEL_CLS_NAME_TAG.format(self._model_cls_name),
            SCVI_VERSION_TAG.format(self._scvi_version),
            ANNDATA_VERSION_TAG.format(self._anndata_version),
        ]
        for m in self._data_modalities:
            tags.append(MODALITY_TAG.format(m))
        for t in self._tissues:
            tags.append(TISSUE_TAG.format(t))
        if self._data_is_annotated is not None:
            tags.append(ANNOTATED_TAG.format(self._data_is_annotated))

        # define the card data, which is the header
        card_data = ModelCardData(
            license=self._license_info,
            library_name=HF_LIBRARY_NAME,
            tags=tags,
        )

        # create the content from the template
        content = template.format(
            card_data=card_data.to_yaml(),
            description=self._description,
            cell_count=self._data_cell_count,
            gene_count=self._data_gene_count,
            model_init_params=json.dumps(self._model_init_params, indent=4),
            model_setup_anndata_args=json.dumps(
                self._model_setup_anndata_args, indent=4
            ),
            model_parent_module=self._model_parent_module,
            large_data_url=self._large_data_url or DEFAULT_NA_FIELD,
            references=self._references,
        )

        # finally create and return the actual card
        # TODO run card.validate()? is it slow?
        return ModelCard(content)

    @property
    def model_card(self) -> ModelCard:
        """Placeholder docstring. TODO complete"""
        return self._model_card
