import pytest
from anndata import AnnData

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager, fields, synthetic_iid
from scvi.data._utils import _get_adata_minify_type
from scvi.model.base import BaseModelClass


class TestModelClass(BaseModelClass):
    @classmethod
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        labels_key: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        **kwargs,
    ):
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            fields.LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            fields.CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            fields.CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
            fields.NumericalObsField(
                REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False
            ),
            fields.CategoricalJointObsField(
                REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys
            ),
            fields.NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]
        # register new fields if the adata is minified
        adata_minify_type = _get_adata_minify_type(adata)
        if adata_minify_type is not None:
            anndata_fields += cls._get_fields_for_adata_minification(adata_minify_type)
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    def train(self):
        pass


def test_deregister_manager():
    adata = synthetic_iid()
    bdata = synthetic_iid()

    # default deregister
    TestModelClass.setup_anndata(adata)
    TestModelClass.setup_anndata(bdata)
    model = TestModelClass(adata)
    adata_manager = model._get_most_recent_anndata_manager(adata)
    bdata_manager = model._get_most_recent_anndata_manager(bdata)

    model.deregister_manager()
    instance_manager_store = TestModelClass._per_instance_manager_store[model.id]
    class_manager_store = model._setup_adata_manager_store
    assert adata_manager.adata_uuid in instance_manager_store
    assert adata_manager.adata_uuid in class_manager_store
    assert bdata_manager.adata_uuid not in class_manager_store

    with pytest.raises(ValueError):
        model.deregister_manager(bdata)

    # deregister with argument
    TestModelClass.setup_anndata(adata)
    TestModelClass.setup_anndata(bdata)
    model = TestModelClass(adata)
    adata_manager = model._get_most_recent_anndata_manager(adata)
    bdata_manager = model._get_most_recent_anndata_manager(bdata)

    model.deregister_manager(adata)
    instance_manager_store = TestModelClass._per_instance_manager_store[model.id]
    class_manager_store = model._setup_adata_manager_store
    assert adata_manager.adata_uuid not in instance_manager_store
    assert adata_manager.adata_uuid not in class_manager_store
    assert bdata_manager.adata_uuid in class_manager_store
    model.deregister_manager(bdata)
    assert bdata_manager.adata_uuid not in class_manager_store

    with pytest.raises(ValueError):
        model.deregister_manager(adata)
    with pytest.raises(ValueError):
        model.deregister_manager(bdata)
