from scvi.external.wscvi import WVAE, WSCVI
from scvi.data import synthetic_iid


def test_wscvi():
    adata = synthetic_iid()
    model = WSCVI(adata=adata, n_latent=5, loss_type="IWELBO", n_particles=25)
    model.train(
        max_epochs=5,
    )

    # Checking that function output shapes make sense
    idx = adata.obs.labels.values == "label_0"
    n_cells = idx.sum()
    scdl = model._make_data_loader(
        adata=adata, indices=idx, batch_size=128
    )
    outs = model._inference_loop(scdl, n_samples=25, transform_batch=None)
    assert outs["log_px_zs"].shape == outs["log_qz"].shape
    assert outs["log_px_zs"].shape == (25 * n_cells, n_cells)

    # Overall 
    outs = model.get_population_scales(
        indices=idx
    )