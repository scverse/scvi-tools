from typing import List, Optional

import torch
from anndata import AnnData
from mudata import MuData

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data._constants import _MODEL_NAME_KEY, _SETUP_ARGS_KEY
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LabelsWithUnlabeledObsField,
    LayerField,
    MuDataCategoricalJointObsField,
    MuDataCategoricalObsField,
    MuDataLayerField,
    MuDataNumericalJointObsField,
    MuDataProteinLayerField,
    NumericalJointObsField,
    ProteinObsmField,
)
from scvi.model import SCVI

use_gpu = torch.cuda.is_available()


def unsupervised_training_one_epoch(
    adata: AnnData,
    run_setup_anndata: bool = True,
    batch_key: Optional[str] = None,
    labels_key: Optional[str] = None,
):
    if run_setup_anndata:
        SCVI.setup_anndata(adata, batch_key=batch_key, labels_key=labels_key)
    m = SCVI(adata)
    m.train(1, train_size=0.4, use_gpu=use_gpu)


def generic_setup_adata_manager(
    adata: AnnData,
    batch_key: Optional[str] = None,
    labels_key: Optional[str] = None,
    categorical_covariate_keys: Optional[List[str]] = None,
    continuous_covariate_keys: Optional[List[str]] = None,
    layer: Optional[str] = None,
    protein_expression_obsm_key: Optional[str] = None,
    protein_names_uns_key: Optional[str] = None,
) -> AnnDataManager:
    setup_args = locals()
    setup_args.pop("adata")
    setup_method_args = {_MODEL_NAME_KEY: "TestModel", _SETUP_ARGS_KEY: setup_args}

    batch_field = CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key)
    anndata_fields = [
        batch_field,
        LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
        CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
        CategoricalJointObsField(
            REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys
        ),
        NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
    ]
    if protein_expression_obsm_key is not None:
        anndata_fields.append(
            ProteinObsmField(
                REGISTRY_KEYS.PROTEIN_EXP_KEY,
                protein_expression_obsm_key,
                use_batch_mask=True,
                batch_field=batch_field,
                colnames_uns_key=protein_names_uns_key,
                is_count_data=True,
            )
        )
    adata_manager = AnnDataManager(
        fields=anndata_fields, setup_method_args=setup_method_args
    )
    adata_manager.register_fields(adata)
    return adata_manager


def scanvi_setup_adata_manager(
    adata: AnnData,
    unlabeled_category: str,
    batch_key: Optional[str] = None,
    labels_key: Optional[str] = None,
) -> AnnDataManager:
    setup_args = locals()
    setup_args.pop("adata")
    setup_method_args = {_MODEL_NAME_KEY: "TestModel", _SETUP_ARGS_KEY: setup_args}
    anndata_fields = [
        CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
        LabelsWithUnlabeledObsField(
            REGISTRY_KEYS.LABELS_KEY, labels_key, unlabeled_category
        ),
    ]
    adata_manager = AnnDataManager(
        fields=anndata_fields, setup_method_args=setup_method_args
    )
    adata_manager.register_fields(adata)
    return adata_manager


def generic_setup_mudata_manager(
    mdata: MuData,
    layer_mod,
    layer: Optional[str] = None,
    batch_mod: Optional[str] = None,
    batch_key: Optional[str] = None,
    categorical_covariate_mod: Optional[str] = None,
    categorical_covariate_keys: Optional[List[str]] = None,
    continuous_covariate_mod: Optional[str] = None,
    continuous_covariate_keys: Optional[List[str]] = None,
    protein_expression_mod: Optional[str] = None,
    protein_expression_layer: Optional[str] = None,
) -> AnnDataManager:
    setup_args = locals()
    setup_args.pop("mdata")
    setup_method_args = {_MODEL_NAME_KEY: "TestModel", _SETUP_ARGS_KEY: setup_args}

    batch_field = MuDataCategoricalObsField(
        REGISTRY_KEYS.BATCH_KEY, batch_key, mod_key=batch_mod
    )
    anndata_fields = [
        MuDataLayerField(
            REGISTRY_KEYS.X_KEY,
            layer,
            mod_key=layer_mod,
            is_count_data=True,
            mod_required=True,
        ),
        batch_field,
        MuDataCategoricalJointObsField(
            REGISTRY_KEYS.CAT_COVS_KEY,
            categorical_covariate_keys,
            mod_key=categorical_covariate_mod,
        ),
        MuDataNumericalJointObsField(
            REGISTRY_KEYS.CONT_COVS_KEY,
            continuous_covariate_keys,
            mod_key=continuous_covariate_mod,
        ),
    ]
    if protein_expression_mod is not None:
        anndata_fields.append(
            MuDataProteinLayerField(
                REGISTRY_KEYS.PROTEIN_EXP_KEY,
                protein_expression_layer,
                mod_key=protein_expression_mod,
                mod_required=True,
                use_batch_mask=True,
                batch_field=batch_field,
            )
        )
    mdata_manager = AnnDataManager(
        fields=anndata_fields, setup_method_args=setup_method_args
    )
    mdata_manager.register_fields(mdata)
    return mdata_manager
