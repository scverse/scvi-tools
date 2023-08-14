import scvi


def test_scanvi_predict_use_posterior_mean():
    adata = scvi.data.synthetic_iid()
    scvi.model.SCANVI.setup_anndata(
        adata, labels_key="labels", unlabeled_category="label_0"
    )

    model = scvi.model.SCANVI(adata)
    model.train(max_epochs=1)

    _ = model.predict(use_posterior_mean=True)
    _ = model.predict(use_posterior_mean=False)
