import inspect
import json
import logging
import os
from typing import List, Optional

import anndata
import torch

logger = logging.getLogger(__name__)


# TODO use pydantic
class HubMetadata:
    """Placeholder docstring. TODO complete"""

    def __init__(
        self,
        tissues: List[str],
        data_modalities: List[str],
        data_cell_count: int,
        data_gene_count: int,
        data_is_annotated: bool,
        model_cls_name: str,
        model_init_params: str,
        model_setup_anndata_args: str,
        scvi_version: str,
        anndata_version: str,
        data_url: Optional[str] = None,
    ):
        self._tissues = tissues
        self._data_modalities = data_modalities
        self._data_cell_count = data_cell_count
        self._data_gene_count = data_gene_count
        self._data_is_annotated = data_is_annotated
        self._data_url = data_url

        self._model_cls_name = model_cls_name
        self._model_init_params = model_init_params
        self._model_setup_anndata_args = model_setup_anndata_args
        self._scvi_version = scvi_version
        self._anndata_version = anndata_version
        # TODO add model criticism metrics under "evaluation metrics" on hugging face

    @classmethod
    def from_dir(
        cls,
        local_dir: str,
        tissues: List[str],
        data_modalities: List[str],
        data_is_annotated: bool,
        anndata_version: str,
    ):
        """Placeholder docstring. TODO complete"""
        adata = anndata.read_h5ad(f"{local_dir}/adata.h5ad")
        data_cell_count = adata.n_obs
        data_gene_count = adata.n_vars
        data_has_full_counts = os.path.isfile(f"{local_dir}/adata_full.h5ad")

        torch_model = torch.load(f"{local_dir}/model.pt")
        attr_dict = torch_model["attr_dict"]
        model_init_params = attr_dict["init_params_"]
        model_cls_name = attr_dict["registry_"]["model_name"]
        model_setup_anndata_args = attr_dict["registry_"]["setup_args"]
        scvi_version = attr_dict["registry_"]["scvi_version"]

        return cls(
            tissues,
            data_modalities,
            data_cell_count,
            data_gene_count,
            data_is_annotated,
            data_has_full_counts,
            model_cls_name,
            model_init_params,
            model_setup_anndata_args,
            scvi_version,
            anndata_version,
        )

    def __repr__(self):
        return f"PTE metadata object with the following properties:\n{self.get_formatted_attrs()}"

    def get_formatted_attrs(self):
        """Placeholder docstring. TODO complete"""
        attributes = inspect.getmembers(self, lambda a: not (inspect.isroutine(a)))

        def skip_attr(a):
            return (
                a[0].startswith("_abc_")
                or a[0].startswith("_")
                or (a[0].startswith("__") and a[0].endswith("__"))
            )

        attributes = [a for a in attributes if not skip_attr(a)]
        attributes = {a[0]: a[1] for a in attributes}
        return json.dumps(attributes, indent=4)

    @property
    def tissues(self) -> List[str]:
        """Placeholder docstring. TODO complete"""
        return self._tissues

    @property
    def data_modalities(self) -> List[str]:
        """Placeholder docstring. TODO complete"""
        return self._data_modalities

    @property
    def data_cell_count(self) -> int:
        """Placeholder docstring. TODO complete"""
        return self._data_cell_count

    @property
    def data_gene_count(self) -> int:
        """Placeholder docstring. TODO complete"""
        return self._data_gene_count

    @property
    def data_is_annotated(self) -> bool:
        """Placeholder docstring. TODO complete"""
        return self._data_is_annotated

    @property
    def data_url(self) -> Optional[str]:
        """Placeholder docstring. TODO complete"""
        return self._data_url

    @property
    def model_cls_name(self) -> str:
        """Placeholder docstring. TODO complete"""
        return self._model_cls_name

    @property
    def model_init_params(self) -> str:
        """Placeholder docstring. TODO complete"""
        return self._model_init_params

    @property
    def model_setup_anndata_args(self) -> str:
        """Placeholder docstring. TODO complete"""
        return self._model_setup_anndata_args

    @property
    def scvi_version(self) -> str:
        """Placeholder docstring. TODO complete"""
        return self._scvi_version

    @property
    def anndata_version(self) -> str:
        """Placeholder docstring. TODO complete"""
        return self._anndata_version
