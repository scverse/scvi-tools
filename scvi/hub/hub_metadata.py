import json
import logging
import os
from pathlib import Path
from typing import List, Optional, Union

import anndata
import rich
import torch
import yaml
from huggingface_hub import ModelCard, ModelCardData
from parse import parse
from rich.markdown import Markdown

from .model_card_template import template

logger = logging.getLogger(__name__)

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


class HubMetadata:
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
        if os.path.isfile(f"{local_dir}/adata.h5ad"):
            adata = anndata.read_h5ad(f"{local_dir}/adata.h5ad", backed=True)
            cell_count = adata.n_obs
            gene_count = adata.n_vars
        else:
            if data_cell_count is None or data_gene_count is None:
                raise ValueError(
                    "No data found on disk. Please provide `data_cell_count` and `data_gene_count`"
                )
            else:
                cell_count = data_cell_count
                gene_count = data_gene_count

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

    @classmethod
    def from_model_card(cls, model_card: Union[ModelCard, str]):
        """Placeholder docstring. TODO complete"""
        # get the content
        if isinstance(model_card, ModelCard):
            content = model_card.content
        elif isinstance(model_card, str):
            content = Path(model_card).read_text()
        else:
            raise TypeError("Unexpected data type for `model_card`")

        # parse the content based on the template
        parser = parse(template, content)
        description = parser["description"]
        data_cell_count = parser["cell_count"]
        data_gene_count = parser["gene_count"]
        model_init_params = json.loads(parser["model_init_params"])
        model_setup_anndata_args = json.loads(parser["model_setup_anndata_args"])
        model_parent_module = parser["model_parent_module"]
        large_data_url = (
            parser["large_data_url"]
            if parser["large_data_url"] != DEFAULT_NA_FIELD
            else None
        )
        references = parser["references"]

        # parse the card_data using yaml
        card_data = yaml.safe_load(parser["card_data"])
        # will look like something like this:
        # {'license': 'cc-by-4.0',
        # 'library_name': 'scvi-tools',
        # 'tags': ['model_cls_name:SCVI',
        # 'scvi_version:0.17.4',
        # 'anndata_version:0.8.0',
        # 'modality:rna',
        # 'modality:atac-seq',
        # 'tissue:eye',
        # 'tissue:spleen']}
        license_info = card_data["license"]
        tags = card_data["tags"]
        model_cls_name = [
            parse(MODEL_CLS_NAME_TAG, t)[0]
            for t in tags
            if parse(MODEL_CLS_NAME_TAG, t)
        ][0]
        scvi_version = [
            parse(SCVI_VERSION_TAG, t)[0] for t in tags if parse(SCVI_VERSION_TAG, t)
        ][0]
        anndata_version = [
            parse(ANNDATA_VERSION_TAG, t)[0]
            for t in tags
            if parse(ANNDATA_VERSION_TAG, t)
        ][0]
        data_modalities = [
            parse(MODALITY_TAG, t)[0] for t in tags if parse(MODALITY_TAG, t)
        ]
        tissues = [parse(TISSUE_TAG, t)[0] for t in tags if parse(TISSUE_TAG, t)]
        data_is_annotated = None
        annotated = [
            parse(ANNOTATED_TAG, t)[0] for t in tags if parse(ANNOTATED_TAG, t)
        ]
        if len(annotated) > 0:
            data_is_annotated = annotated[0]

        # instantiate the class
        return cls(
            license_info,
            data_cell_count,
            data_gene_count,
            model_cls_name,
            model_init_params,
            model_setup_anndata_args,
            scvi_version,
            anndata_version,
            data_modalities,
            tissues,
            data_is_annotated,
            large_data_url,
            model_parent_module,
            description,
            references,
        )

    def __repr__(self):
        print("HubMetadata wrapping the following ModelCard 👇")
        rich.print(Markdown(self.model_card.content))
        return ""

    @property
    def model_card(self) -> ModelCard:
        """Placeholder docstring. TODO complete"""
        return self._model_card

    @property
    def large_data_url(self) -> Optional[str]:
        """Placeholder docstring. TODO complete"""
        return self._large_data_url

    @property
    def model_parent_module(self) -> Optional[str]:
        """Placeholder docstring. TODO complete"""
        return self._model_parent_module
