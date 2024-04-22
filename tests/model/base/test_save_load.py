import numpy as np
from pytest import raises

from scvi.data import synthetic_iid
from scvi.model.base._save_load import _get_var_names


def test_get_var_names_anndata(n_genes: int = 50):
    adata = synthetic_iid(n_genes=n_genes)

    var_names = _get_var_names(adata)
    assert isinstance(var_names, np.ndarray)
    assert len(var_names) == n_genes


def test_get_var_names_mudata(n_genes: int = 50, n_proteins: int = 10):
    mdata = synthetic_iid(n_genes=n_genes, n_proteins=n_proteins, n_regions=0, return_mudata=True)

    var_names = _get_var_names(mdata)
    assert isinstance(var_names, dict)
    assert all(mod_key in mdata.mod for mod_key in var_names)
    assert all(isinstance(v, np.ndarray) for v in var_names.values())
    assert len(var_names["rna"]) == n_genes
    assert len(var_names["protein_expression"]) == n_proteins

    legacy_var_names = _get_var_names(mdata, legacy_mudata_format=True)
    assert isinstance(legacy_var_names, np.ndarray)
    assert len(legacy_var_names) == n_genes + n_proteins


def test_get_var_names_type():
    with raises(TypeError):
        _get_var_names({})
