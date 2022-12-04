import json
import logging
from pathlib import Path
from typing import List, Optional

import anndata
import torch
from huggingface_hub import ModelCard, ModelCardData

logger = logging.getLogger(__name__)

HF_LIBRARY_NAME = "scvi-tools"
MODEL_CARD_TEMPLATE_FILE = "model_card_template.md"


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
        data_url: Optional[str] = None,
        description: str = "",
    ):
        self._data_cell_count = data_cell_count
        self._data_gene_count = data_gene_count
        self._data_modalities = data_modalities or []
        self._tissues = tissues or []
        self._data_is_annotated = data_is_annotated
        self._data_url = data_url

        self._model_cls_name = model_cls_name
        self._model_init_params = model_init_params
        self._model_setup_anndata_args = model_setup_anndata_args

        self._license_info = license_info
        self._scvi_version = scvi_version
        self._anndata_version = anndata_version
        self._description = description

        # TODO add model criticism metrics under "evaluation metrics" on hugging face

        self._model_card = self._to_model_card()

    @classmethod
    def from_dir(
        cls,
        local_dir: str,
        license_info: str,
        anndata_version: str,
        tissues: Optional[List[str]] = None,
        data_modalities: Optional[List[str]] = None,
        data_is_annotated: Optional[bool] = None,
        description: str = "",
    ):
        """Placeholder docstring. TODO complete"""
        adata = anndata.read_h5ad(f"{local_dir}/adata.h5ad")
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
            data_modalities,
            tissues,
            data_is_annotated,
            None,
            description,
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
            tags.append(f"modality:{m}")
        for t in self._tissues:
            tags.append(f"tissue:{t}")
        if self._data_is_annotated is not None:
            tags.append(f"annotated:{self._data_is_annotated  }")

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
        )

        # finally create and return the actual card
        return ModelCard(content)

    def __repr__(self):
        return f"HubMetadata wrapping the following ModelCard:\n{self.model_card}"

    @property
    def model_card(self) -> ModelCard:
        """Placeholder docstring. TODO complete"""
        return self._model_card
