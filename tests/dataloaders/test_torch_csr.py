import scvi


def test_torch_csr_scvi():
    adata = scvi.data.synthetic_iid()
    scvi.model.SCVI.setup_anndata(adata)

    model = scvi.model.SCVI(adata)
    model.train(max_epochs=1, transfer_torch_csr=True)
