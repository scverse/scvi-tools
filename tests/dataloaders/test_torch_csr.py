import scvi


def test_torch_csr_scvi():
    adata = scvi.data.synthetic_iid(sparse_format="csr_array")
    scvi.model.SCVI.setup_anndata(adata)

    model = scvi.model.SCVI(adata)
    model.train(max_epochs=1, load_sparse_tensor=True)
