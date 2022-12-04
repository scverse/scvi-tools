import json
import logging
from pathlib import Path
from typing import List, Optional, Union

import anndata
import torch
import yaml
from huggingface_hub import ModelCard, ModelCardData
from parse import parse, search

logger = logging.getLogger(__name__)

HF_LIBRARY_NAME = "scvi-tools"
MODEL_CARD_TEMPLATE_FILE = "model_card_template.md"
DEFAULT_MISSING_FIELD = "To be added..."
DEFAULT_NA_FIELD = "N/A"
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
        model_init_params: str,
        model_setup_anndata_args: str,
        scvi_version: str,
        anndata_version: str,
        data_modalities: Optional[List[str]] = None,
        tissues: Optional[List[str]] = None,
        data_is_annotated: Optional[bool] = None,
        large_data_url: Optional[str] = None,
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

        self._license_info = license_info
        self._scvi_version = scvi_version
        self._anndata_version = anndata_version
        self._description = description
        self._references = references

        # TODO add model criticism metrics under "evaluation metrics" on hugging face

        self._model_card = self._to_model_card()

    @classmethod
    def from_dir(
        cls,
        local_dir: str,
        license_info: str,
        anndata_version: str,
        **kwargs,
    ):
        """Placeholder docstring. TODO complete"""
        adata = anndata.read_h5ad(f"{local_dir}/adata.h5ad", backed=True)
        data_cell_count = adata.n_obs
        data_gene_count = adata.n_vars

        torch_model = torch.load(f"{local_dir}/model.pt")
        attr_dict = torch_model["attr_dict"]
        model_init_params = attr_dict["init_params_"]
        model_cls_name = attr_dict["registry_"]["model_name"]
        model_setup_anndata_args = attr_dict["registry_"]["setup_args"]
        scvi_version = attr_dict["registry_"]["scvi_version"]

        return cls(
            license_info,
            data_cell_count,
            data_gene_count,
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
            f"model_cls_name:{self._model_cls_name}",
            f"scvi_version:{self._scvi_version}",
            f"anndata_version:{self._anndata_version}",
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
            library_name="scvi-tools",
            tags=tags,
        )

        # create the content from the template
        template = (Path(__file__).parent / MODEL_CARD_TEMPLATE_FILE).read_text()
        content = template.format(
            card_data=card_data.to_yaml(),
            description=self._description,
            cell_count=self._data_cell_count,
            gene_count=self._data_gene_count,
            model_init_params=json.dumps(self._model_init_params, indent=4),
            model_setup_anndata_args=json.dumps(
                self._model_setup_anndata_args, indent=4
            ),
            large_data_url=self._large_data_url or DEFAULT_NA_FIELD,
            references=self._references,
        )

        # finally create and return the actual card
        return ModelCard(content)

    @classmethod
    def from_model_card(cls, model_card: Union[ModelCard, str]):
        """Placeholder docstring. TODO complete"""
        # get the template
        template = (Path(__file__).parent / MODEL_CARD_TEMPLATE_FILE).read_text()

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
        model_init_params = parser["model_init_params"]
        model_setup_anndata_args = parser["model_setup_anndata_args"]
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
        model_cls_name = tags["model_cls_name"]
        scvi_version = tags["scvi_version"]
        anndata_version = tags["anndata_version"]
        data_modalities = [
            search(MODALITY_TAG, t)[0] for t in tags if search(MODALITY_TAG, t)
        ]
        tissues = [search(TISSUE_TAG, t)[0] for t in tags if search(TISSUE_TAG, t)]
        data_is_annotated = None
        annotated = [
            search(ANNOTATED_TAG, t)[0] for t in tags if search(ANNOTATED_TAG, t)
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
            description,
            references,
        )

    def __repr__(self):
        return f"HubMetadata wrapping the following ModelCard:\n{self.model_card}"

    @property
    def model_card(self) -> ModelCard:
        """Placeholder docstring. TODO complete"""
        return self._model_card

    @property
    def large_data_url(self) -> Optional[str]:
        """Placeholder docstring. TODO complete"""
        return self._large_data_url
