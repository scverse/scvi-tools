from torch.nn import Linear

from scvi.data import synthetic_iid
from scvi.external import POISSONVI


def test_poissonvi():
    adata = synthetic_iid(batch_size=100)
    POISSONVI.setup_anndata(adata, batch_key="batch")
    model = POISSONVI(adata)
    model.train(max_epochs=1)
    model.get_latent_representation()
    model.get_normalized_accessibility()
    model.get_region_factors()
    model.get_normalized_expression()
    model.get_normalized_expression(transform_batch="batch_1")
    model.get_normalized_expression(n_samples=2)


def test_poissonvi_default_params():
    from scvi.model import PEAKVI

    adata = synthetic_iid(batch_size=100)
    POISSONVI.setup_anndata(adata)
    PEAKVI.setup_anndata(adata)
    poissonvi = POISSONVI(adata)
    peakvi = PEAKVI(adata)

    assert poissonvi.module.n_latent == peakvi.module.n_latent
    assert poissonvi.module.latent_distribution == peakvi.module.latent_distribution
    poisson_encoder = poissonvi.module.z_encoder.encoder
    poisson_mean_encoder = poissonvi.module.z_encoder.mean_encoder
    poisson_decoder = poissonvi.module.decoder.px_decoder
    assert len(poisson_encoder.fc_layers) == peakvi.module.n_layers_encoder
    assert len(poisson_decoder.fc_layers) == peakvi.module.n_layers_encoder
    assert poisson_encoder.fc_layers[-1][0].in_features == peakvi.module.n_hidden
    assert poisson_decoder.fc_layers[-1][0].in_features == peakvi.module.n_hidden
    assert poisson_mean_encoder.out_features == peakvi.module.n_latent
    assert poisson_decoder.fc_layers[0][0].in_features == peakvi.module.n_latent


def test_poissonvi_non_default_params(
    n_hidden: int = 50,
    n_latent: int = 5,
    n_layers: int = 2,
    dropout_rate: float = 0.4,
    latent_distribution="ln",
):
    adata = synthetic_iid(batch_size=100)
    POISSONVI.setup_anndata(adata)
    model = POISSONVI(
        adata,
        n_hidden=n_hidden,
        n_latent=n_latent,
        n_layers=n_layers,
        dropout_rate=dropout_rate,
        latent_distribution=latent_distribution,
    )

    assert model.module.n_latent == n_latent
    assert model.module.latent_distribution == latent_distribution

    encoder = model.module.z_encoder.encoder
    assert len(encoder.fc_layers) == n_layers
    linear = encoder.fc_layers[-1][0]
    assert isinstance(linear, Linear)
    assert linear.in_features == n_hidden
    mean_encoder = model.module.z_encoder.mean_encoder
    assert isinstance(mean_encoder, Linear)
    assert mean_encoder.out_features == n_latent

    model.train(max_epochs=1)
    assert model.get_latent_representation().shape[1] == n_latent
