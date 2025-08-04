from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data._utils import _get_adata_minify_type
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LabelsWithUnlabeledObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
    ObsmField,
)
from scvi.model._scanvi import SCANVI
from scvi.utils import setup_anndata_dsp

from ._module import HSCANVAE
from ._utils import process_adata_ontology

NODE_KEY = "hscanvi_co_node"
NODE_IDX_KEY = "hscanvi_co_node_idx"
MULTILABEL_KEY = "hscanvi_co_multilabel_onehot"

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from anndata import AnnData


class HSCANVI(SCANVI):
    _module_cls = HSCANVAE

    def __init__(
        self,
        adata: AnnData,
        **model_kwargs,
    ):
        super().__init__(adata, **model_kwargs)
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        labels_key: str,
        unlabeled_category: str,
        layer: str | None = None,
        batch_key: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        use_minified: bool = True,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_labels_key)s
        %(param_unlabeled_category)s
        %(param_layer)s
        %(param_batch_key)s
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        use_minified
            If True, will register the minified version of the adata if possible.
        """
        process_adata_ontology(
            adata,
            cell_type_key=labels_key,
            added_node_key=NODE_KEY,
            added_node_idx_key=NODE_IDX_KEY,
            added_multilabel_key=MULTILABEL_KEY,
        )
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            LabelsWithUnlabeledObsField(REGISTRY_KEYS.LABELS_KEY, labels_key, unlabeled_category),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
            ObsmField("multilabels", MULTILABEL_KEY),
        ]
        # register new fields if the adata is minified
        if adata:
            adata_minify_type = _get_adata_minify_type(adata)
            if adata_minify_type is not None and use_minified:
                anndata_fields += cls._get_fields_for_adata_minification(adata_minify_type)
            adata_manager = AnnDataManager(
                fields=anndata_fields, setup_method_args=setup_method_args
            )
            adata_manager.register_fields(adata, **kwargs)
            cls.register_manager(adata_manager)
