import numpy as np
import pandas as pd
import pytest
import torch
from anndata import AnnData
from mudata import MuData

from scvi.external.diagvi import DIAGVI


# =============================================================================
# Constants and fixtures
# =============================================================================

# Constants for test data dimensions
N_OBS_SEQ = 100
N_OBS_SPATIAL = 80
N_VARS = 50
N_LABELS = 5  # Number of cell types for semi-supervised tests


@pytest.fixture(scope="module")
def adata_seq():
    X = np.random.poisson(1.0, size=(N_OBS_SEQ, N_VARS))
    obs = {"batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=N_OBS_SEQ))}
    var = pd.DataFrame(index=[f"gene{i}" for i in range(N_VARS)])
    return AnnData(X=X, obs=obs, var=var)


@pytest.fixture(scope="module")
def adata_spatial():
    X = np.random.poisson(1.0, size=(N_OBS_SPATIAL, N_VARS))
    obs = {"batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=N_OBS_SPATIAL))}
    var = pd.DataFrame(index=[f"gene{i}" for i in range(N_VARS)])
    return AnnData(X=X, obs=obs, var=var)


@pytest.fixture(scope="module")
def trained_model(adata_seq, adata_spatial):
    """Module-scoped trained model to avoid redundant training across tests."""
    DIAGVI.setup_anndata(adata_seq, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nb")
    model = DIAGVI({"diss": adata_seq, "spatial": adata_spatial})
    model.train(max_epochs=1, batch_size=16)
    return model, adata_seq, adata_spatial


def make_model(adata_seq, adata_spatial):
    DIAGVI.setup_anndata(adata_seq, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nb")
    return DIAGVI({"diss": adata_seq, "spatial": adata_spatial})


# =============================================================================
# Basic model functionality tests
# =============================================================================


def test_diagvi_initialization(adata_seq, adata_spatial):
    model = make_model(adata_seq, adata_spatial)
    assert model.module is not None
    assert "diss" in model.adatas
    assert "spatial" in model.adatas


def test_diagvi_training(adata_seq, adata_spatial):
    model = make_model(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)
    assert model.is_trained_ is True


def test_diagvi_latent_representation(adata_seq, adata_spatial):
    model = make_model(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)
    latents = model.get_latent_representation()
    assert latents["diss"].shape[0] == adata_seq.shape[0]
    assert latents["spatial"].shape[0] == adata_spatial.shape[0]


def test_diagvi_generative_model(adata_seq, adata_spatial):
    model = make_model(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)
    tensors = {
        "X": torch.tensor(adata_seq.X, dtype=torch.float32),
        "batch": torch.tensor(adata_seq.obs["batch"].cat.codes.values, dtype=torch.long),
    }
    inference_outputs = model.module.inference(tensors["X"], mode="diss")
    generative_outputs = model.module.generative(
        inference_outputs["z"],
        inference_outputs["library"],
        batch_index=tensors["batch"],
        v=inference_outputs["v"],
        mode="diss",
    )
    assert "px" in generative_outputs or "px_rate" in generative_outputs
    assert "pz" in generative_outputs


def test_diagvi_loss(adata_seq, adata_spatial):
    model = make_model(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)
    device = next(model.module.parameters()).device
    tensors = {
        "X": torch.tensor(adata_seq.X, dtype=torch.float32, device=device),
        "batch": torch.tensor(
            adata_seq.obs["batch"].cat.codes.values, dtype=torch.long, device=device
        ),
    }
    inference_outputs = model.module.inference(tensors["X"], mode="diss")
    generative_outputs = model.module.generative(
        inference_outputs["z"],
        inference_outputs["library"],
        batch_index=tensors["batch"],
        v=inference_outputs["v"],
        mode="diss",
    )
    loss = model.module.loss(
        tensors,
        inference_outputs,
        generative_outputs,
        lam_kl=1.0, lam_data=0.1,
        mode="diss")
    assert loss.loss is not None
    assert loss.reconstruction_loss is not None
    assert loss.kl_local is not None


def test_diagvi_get_imputed_values(adata_seq, adata_spatial):
    model = make_model(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)
    # Impute spatial from diss
    imputed = model.get_imputed_values(query_name="spatial")
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape[0] == adata_spatial.shape[0]
    assert imputed.shape[1] == adata_seq.shape[1]
    # Impute diss from spatial
    imputed2 = model.get_imputed_values(query_name="diss")
    assert isinstance(imputed2, np.ndarray)
    assert imputed2.shape[0] == adata_seq.shape[0]
    assert imputed2.shape[1] == adata_spatial.shape[1]


def test_diagvi_save_and_load(adata_seq, adata_spatial, tmp_path):
    model = make_model(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)
    save_path = tmp_path / "diagvi_test_model.pt"
    model.save(save_path, save_anndata=True)

    # Load the model
    loaded_model = DIAGVI.load(save_path)

    # compare state dicts
    model_state = model.module.state_dict()
    loaded_state = loaded_model.module.state_dict()
    for key in model_state:
        assert key in loaded_state
        assert torch.allclose(model_state[key], loaded_state[key], atol=1e-6), f"Mismatch in {key}"

    # compare output equivalence
    latents = loaded_model.get_latent_representation()
    assert latents["diss"].shape[0] == adata_seq.shape[0]
    assert latents["spatial"].shape[0] == adata_spatial.shape[0]


# =============================================================================
# Tests for get_imputed_values
# =============================================================================


@pytest.mark.parametrize("reference_batch", [0, 1, "batch1", "batch2"])
def test_get_imputed_values_reference_batch_scalar(trained_model, reference_batch):
    """Test get_imputed_values with scalar reference_batch (int or string)."""
    model, adata_seq, adata_spatial = trained_model

    # Impute from spatial to diss
    imputed = model.get_imputed_values(query_name="spatial", reference_batch=reference_batch)
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape == (N_OBS_SPATIAL, N_VARS)

    # Impute from diss to spatial
    imputed2 = model.get_imputed_values(query_name="diss", reference_batch=reference_batch)
    assert isinstance(imputed2, np.ndarray)
    assert imputed2.shape == (N_OBS_SEQ, N_VARS)


@pytest.mark.parametrize(
    "batch_array_type",
    ["int_array", "string_array"],
)
def test_get_imputed_values_reference_batch_array(trained_model, batch_array_type):
    """Test get_imputed_values with array-like reference_batch."""
    model, adata_seq, adata_spatial = trained_model

    if batch_array_type == "int_array":
        # Integer array matching spatial obs count
        reference_batch = np.random.choice([0, 1], size=N_OBS_SPATIAL)
    else:
        # String array matching spatial obs count
        reference_batch = np.random.choice(["batch1", "batch2"], size=N_OBS_SPATIAL)

    imputed = model.get_imputed_values(query_name="spatial", reference_batch=reference_batch)
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape == (N_OBS_SPATIAL, N_VARS)


@pytest.mark.parametrize("reference_libsize", [100.0, 1000.0, 10000.0])
def test_get_imputed_values_reference_libsize_scalar(trained_model, reference_libsize):
    """Test get_imputed_values with scalar reference_libsize."""
    model, adata_seq, adata_spatial = trained_model

    imputed = model.get_imputed_values(query_name="spatial", reference_libsize=reference_libsize)
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape == (N_OBS_SPATIAL, N_VARS)


def test_get_imputed_values_reference_libsize_array(trained_model):
    """Test get_imputed_values with array reference_libsize."""
    model, adata_seq, adata_spatial = trained_model

    # Array of library sizes matching obs count
    reference_libsize = np.random.uniform(500, 2000, size=N_OBS_SPATIAL)

    imputed = model.get_imputed_values(query_name="spatial", reference_libsize=reference_libsize)
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape == (N_OBS_SPATIAL, N_VARS)


def test_get_imputed_values_combined_batch_and_libsize(trained_model):
    """Test get_imputed_values with both reference_batch and reference_libsize specified."""
    model, adata_seq, adata_spatial = trained_model

    # Use scalar values
    imputed = model.get_imputed_values(
        query_name="spatial", reference_batch=1, reference_libsize=1000.0
    )
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape == (N_OBS_SPATIAL, N_VARS)

    # Use array values
    reference_batch_arr = np.random.choice([0, 1], size=N_OBS_SPATIAL)
    reference_libsize_arr = np.random.uniform(500, 2000, size=N_OBS_SPATIAL)

    imputed2 = model.get_imputed_values(
        query_name="spatial",
        reference_batch=reference_batch_arr,
        reference_libsize=reference_libsize_arr,
    )
    assert isinstance(imputed2, np.ndarray)
    assert imputed2.shape == (N_OBS_SPATIAL, N_VARS)

    # Mixed: scalar batch, array libsize
    imputed3 = model.get_imputed_values(
        query_name="spatial",
        reference_batch="batch1",
        reference_libsize=reference_libsize_arr,
    )
    assert isinstance(imputed3, np.ndarray)
    assert imputed3.shape == (N_OBS_SPATIAL, N_VARS)


@pytest.mark.parametrize(
    ("error_case", "error_kwargs"),
    [
        ("wrong_size_batch", {"reference_batch": np.array([0, 1, 0])}),  # Wrong size
        ("wrong_size_libsize", {"reference_libsize": np.array([100.0, 200.0, 300.0])}),  # Wrong size
        ("2d_libsize", {"reference_libsize": np.array([[100.0], [200.0], [300.0]])}),  # 2D array
    ],
)
def test_get_imputed_values_errors(trained_model, error_case, error_kwargs):
    """Test that get_imputed_values raises ValueError for invalid inputs."""
    model, adata_seq, adata_spatial = trained_model

    with pytest.raises(ValueError):
        model.get_imputed_values(query_name="spatial", **error_kwargs)


# =============================================================================
# Tests for MuData support
# =============================================================================


@pytest.fixture(scope="module")
def mudata_fixture():
    """Create a MuData object with two modalities for testing."""
    # Create first modality (e.g., RNA)
    X_rna = np.random.poisson(1.0, size=(N_OBS_SEQ, N_VARS))
    obs_rna = pd.DataFrame(
        {"batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=N_OBS_SEQ))}
    )
    var_rna = pd.DataFrame(index=[f"gene{i}" for i in range(N_VARS)])
    adata_rna = AnnData(X=X_rna, obs=obs_rna, var=var_rna)

    # Create second modality (e.g., spatial)
    X_spatial = np.random.poisson(1.0, size=(N_OBS_SPATIAL, N_VARS))
    obs_spatial = pd.DataFrame(
        {"batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=N_OBS_SPATIAL))}
    )
    var_spatial = pd.DataFrame(index=[f"gene{i}" for i in range(N_VARS)])
    adata_spatial = AnnData(X=X_spatial, obs=obs_spatial, var=var_spatial)

    # Create MuData
    mdata = MuData({"rna": adata_rna, "spatial": adata_spatial})
    return mdata


def test_diagvi_mudata_setup(mudata_fixture):
    """Test that setup_mudata correctly registers both modalities."""
    mdata = mudata_fixture

    DIAGVI.setup_mudata(
        mdata,
        modalities=["rna", "spatial"],
        batch_key="batch",
        likelihood="nb",
    )

    # Check that diagvi-specific keys are set in uns
    assert "diagvi_likelihood" in mdata.mod["rna"].uns
    assert "diagvi_likelihood" in mdata.mod["spatial"].uns
    assert mdata.mod["rna"].uns["diagvi_likelihood"] == "nb"
    assert mdata.mod["spatial"].uns["diagvi_likelihood"] == "nb"


def test_diagvi_mudata_initialization(mudata_fixture):
    """Test that DIAGVI can be initialized with a MuData object."""
    mdata = mudata_fixture

    DIAGVI.setup_mudata(
        mdata,
        modalities=["rna", "spatial"],
        batch_key="batch",
        likelihood="nb",
    )

    # Initialize model with MuData directly
    model = DIAGVI(mdata)

    assert model.module is not None
    assert "rna" in model.adatas
    assert "spatial" in model.adatas
    assert len(model.input_names) == 2


def test_diagvi_mudata_training(mudata_fixture):
    """Test that DIAGVI can be trained when initialized with MuData."""
    mdata = mudata_fixture

    DIAGVI.setup_mudata(
        mdata,
        modalities=["rna", "spatial"],
        batch_key="batch",
        likelihood="nb",
    )

    model = DIAGVI(mdata)
    model.train(max_epochs=1, batch_size=16)

    assert model.is_trained_ is True


def test_diagvi_mudata_latent_representation(mudata_fixture):
    """Test latent representation extraction with MuData-initialized model."""
    mdata = mudata_fixture

    DIAGVI.setup_mudata(
        mdata,
        modalities=["rna", "spatial"],
        batch_key="batch",
        likelihood="nb",
    )

    model = DIAGVI(mdata)
    model.train(max_epochs=1, batch_size=16)

    latents = model.get_latent_representation()

    assert "rna" in latents
    assert "spatial" in latents
    assert latents["rna"].shape[0] == N_OBS_SEQ
    assert latents["spatial"].shape[0] == N_OBS_SPATIAL


def test_diagvi_mudata_imputation(mudata_fixture):
    """Test imputation with MuData-initialized model."""
    mdata = mudata_fixture

    DIAGVI.setup_mudata(
        mdata,
        modalities=["rna", "spatial"],
        batch_key="batch",
        likelihood="nb",
    )

    model = DIAGVI(mdata)
    model.train(max_epochs=1, batch_size=16)

    # Impute from spatial to rna
    imputed = model.get_imputed_values(query_name="spatial")
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape[0] == N_OBS_SPATIAL

    # Impute from rna to spatial
    imputed2 = model.get_imputed_values(query_name="rna")
    assert isinstance(imputed2, np.ndarray)
    assert imputed2.shape[0] == N_OBS_SEQ


def test_diagvi_mudata_with_different_params_per_modality(mudata_fixture):
    """Test setup_mudata with modality-specific parameters passed as dicts."""
    mdata = mudata_fixture

    DIAGVI.setup_mudata(
        mdata,
        modalities=["rna", "spatial"],
        batch_key={"rna": "batch", "spatial": "batch"},
        likelihood={"rna": "nb", "spatial": "zinb"},
    )

    assert mdata.mod["rna"].uns["diagvi_likelihood"] == "nb"
    assert mdata.mod["spatial"].uns["diagvi_likelihood"] == "zinb"

    model = DIAGVI(mdata)
    assert model.module is not None


# =============================================================================
# GMM prior and semi-supervised fixtures and tests
# =============================================================================


@pytest.fixture(scope="module")
def adata_seq_with_labels():
    """AnnData with cell type labels for semi-supervised testing."""
    X = np.random.poisson(1.0, size=(N_OBS_SEQ, N_VARS))
    labels = np.random.choice([f"celltype_{i}" for i in range(N_LABELS)], size=N_OBS_SEQ)
    obs = pd.DataFrame(
        {
            "batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=N_OBS_SEQ)),
            "cell_type": pd.Categorical(labels),
        }
    )
    var = pd.DataFrame(index=[f"gene{i}" for i in range(N_VARS)])
    return AnnData(X=X, obs=obs, var=var)


@pytest.fixture(scope="module")
def adata_spatial_with_labels():
    """AnnData with cell type labels for semi-supervised testing."""
    X = np.random.poisson(1.0, size=(N_OBS_SPATIAL, N_VARS))
    labels = np.random.choice([f"celltype_{i}" for i in range(N_LABELS)], size=N_OBS_SPATIAL)
    obs = pd.DataFrame(
        {
            "batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=N_OBS_SPATIAL)),
            "cell_type": pd.Categorical(labels),
        }
    )
    var = pd.DataFrame(index=[f"gene{i}" for i in range(N_VARS)])
    return AnnData(X=X, obs=obs, var=var)


def test_diagvi_gmm_prior_custom_components(adata_seq, adata_spatial):
    """Test DIAGVI with custom number of GMM components."""
    n_components = 7
    DIAGVI.setup_anndata(
        adata_seq,
        batch_key="batch",
        likelihood="nb",
        gmm_prior=True,
        n_mixture_components=n_components,
    )
    DIAGVI.setup_anndata(
        adata_spatial,
        batch_key="batch",
        likelihood="nb",
        gmm_prior=True,
        n_mixture_components=n_components,
    )

    model = DIAGVI({"diss": adata_seq, "spatial": adata_spatial})

    assert model.module.gmm_means["diss"].shape[0] == n_components
    assert model.module.gmm_means["spatial"].shape[0] == n_components


def test_diagvi_semi_supervised_one_modality(adata_seq_with_labels, adata_spatial):
    """Test semi-supervised on only one modality."""
    # Semi-supervised is automatically enabled when labels_key is provided
    DIAGVI.setup_anndata(
        adata_seq_with_labels,
        batch_key="batch",
        labels_key="cell_type",
        likelihood="nb",
        gmm_prior=True,
    )
    # No labels_key means no semi-supervised mode
    DIAGVI.setup_anndata(
        adata_spatial, batch_key="batch", likelihood="nb", gmm_prior=False
    )

    model = DIAGVI({"diss": adata_seq_with_labels, "spatial": adata_spatial})

    # First modality should have GMM prior and classifier
    assert model.module.use_gmm_prior["diss"] is True
    assert model.module.classifier_0 is not None
    assert model.module.gmm_means["diss"].shape[0] == model.module.n_labels["diss"]

    # Second modality should have standard normal prior and no classifier
    assert model.module.use_gmm_prior["spatial"] is False
    assert model.module.classifier_1 is None
    assert "spatial" not in model.module.gmm_means

    # Should still train successfully
    model.train(max_epochs=1, batch_size=16)
    assert model.is_trained_ is True


# =============================================================================
# Tests for automatic semi_supervised handling
# =============================================================================


def test_semi_supervised_auto_enabled_with_labels_key(adata_seq_with_labels):
    """Test that semi_supervised is automatically enabled when labels_key is provided."""
    DIAGVI.setup_anndata(
        adata_seq_with_labels,
        batch_key="batch",
        labels_key="cell_type",
        likelihood="nb",
    )
    
    assert adata_seq_with_labels.uns["diagvi_semi_supervised"] is True
    # n_mixture_components should equal number of unique labels
    expected_n_labels = adata_seq_with_labels.obs["cell_type"].nunique()
    assert adata_seq_with_labels.uns["diagvi_n_mixture_components"] == expected_n_labels


def test_semi_supervised_auto_disabled_without_labels_key(adata_seq):
    """Test that semi_supervised is disabled when no labels_key is provided."""
    DIAGVI.setup_anndata(
        adata_seq,
        batch_key="batch",
        likelihood="nb",
        n_mixture_components=15,
    )
    
    assert adata_seq.uns["diagvi_semi_supervised"] is False
    # n_mixture_components should be the user-specified value
    assert adata_seq.uns["diagvi_n_mixture_components"] == 15


def test_n_mixture_components_ignored_with_labels_key_warning(adata_seq_with_labels):
    """Test that warning is raised when n_mixture_components differs from label count."""
    expected_n_labels = adata_seq_with_labels.obs["cell_type"].nunique()
    
    with pytest.warns(UserWarning, match="n_mixture_components=99 is ignored"):
        DIAGVI.setup_anndata(
            adata_seq_with_labels,
            batch_key="batch",
            labels_key="cell_type",
            likelihood="nb",
            n_mixture_components=99,  # This should be ignored
        )
    
    # Should use label count, not user-specified value
    assert adata_seq_with_labels.uns["diagvi_n_mixture_components"] == expected_n_labels


def test_n_mixture_components_no_warning_when_matching(adata_seq_with_labels):
    """Test that no warning when n_mixture_components matches label count."""
    expected_n_labels = adata_seq_with_labels.obs["cell_type"].nunique()
    
    # Should not warn if the values match
    import warnings as warn_module
    with warn_module.catch_warnings():
        warn_module.simplefilter("error")  # Turn warnings into errors
        DIAGVI.setup_anndata(
            adata_seq_with_labels,
            batch_key="batch",
            labels_key="cell_type",
            likelihood="nb",
            n_mixture_components=expected_n_labels,
        )


def test_n_mixture_components_no_warning_when_default(adata_seq_with_labels):
    """Test that no warning when n_mixture_components is default value (10)."""
    import warnings as warn_module
    with warn_module.catch_warnings():
        warn_module.simplefilter("error")  # Turn warnings into errors
        DIAGVI.setup_anndata(
            adata_seq_with_labels,
            batch_key="batch",
            labels_key="cell_type",
            likelihood="nb",
            # n_mixture_components defaults to 10, should not warn
        )


# =============================================================================
# Tests for construct_custom_guidance_graph
# =============================================================================


@pytest.fixture(scope="module")
def guidance_graph_adatas():
    """AnnData objects for testing construct_custom_guidance_graph."""
    n_vars_diss = 20
    n_vars_spatial = 15
    adata_diss = AnnData(
        X=np.random.poisson(1.0, size=(50, n_vars_diss)),
        obs={"batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=50))},
        var=pd.DataFrame(index=[f"gene{i}" for i in range(n_vars_diss)]),
    )
    adata_spatial = AnnData(
        X=np.random.poisson(1.0, size=(40, n_vars_spatial)),
        obs={"batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=40))},
        var=pd.DataFrame(index=[f"gene{i}" for i in range(n_vars_spatial)]),
    )
    return adata_diss, adata_spatial


def test_construct_custom_guidance_graph_basic(guidance_graph_adatas):
    """Test basic construction of custom guidance graph."""
    adata_diss, adata_spatial = guidance_graph_adatas
    DIAGVI.setup_anndata(adata_diss, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nb")

    # Create mapping between some genes
    mapping_df = pd.DataFrame(
        {
            "diss": ["gene0", "gene1", "gene2"],
            "spatial": ["gene0", "gene1", "gene2"],
        }
    )

    graph = DIAGVI.construct_custom_guidance_graph(
        input_dict={"diss": adata_diss, "spatial": adata_spatial},
        mapping_df=mapping_df,
    )

    # Check that graph is a torch_geometric Data object
    assert hasattr(graph, "edge_index")
    assert hasattr(graph, "edge_weight")
    assert hasattr(graph, "edge_sign")
    assert hasattr(graph, "x")
    assert hasattr(graph, "diss_indices")
    assert hasattr(graph, "spatial_indices")

    # Check node count: diss (20) + spatial (15) = 35
    assert graph.x.shape[0] == 20 + 15


def test_construct_custom_guidance_graph_edge_count(guidance_graph_adatas):
    """Test that edge count is correct: bidirectional edges + self-loops."""
    adata_diss, adata_spatial = guidance_graph_adatas
    DIAGVI.setup_anndata(adata_diss, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nb")

    # 3 mapped pairs
    mapping_df = pd.DataFrame(
        {
            "diss": ["gene0", "gene1", "gene2"],
            "spatial": ["gene0", "gene1", "gene2"],
        }
    )

    graph = DIAGVI.construct_custom_guidance_graph(
        input_dict={"diss": adata_diss, "spatial": adata_spatial},
        mapping_df=mapping_df,
    )

    n_nodes = 20 + 15
    n_mapped_pairs = 3
    # Bidirectional edges: 2 * n_mapped_pairs
    # Self-loops: n_nodes
    expected_edges = 2 * n_mapped_pairs + n_nodes

    assert graph.edge_index.shape[1] == expected_edges


def test_construct_custom_guidance_graph_custom_weight_sign(guidance_graph_adatas):
    """Test custom weight and sign parameters."""
    adata_diss, adata_spatial = guidance_graph_adatas
    DIAGVI.setup_anndata(adata_diss, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nb")

    mapping_df = pd.DataFrame(
        {
            "diss": ["gene0", "gene1"],
            "spatial": ["gene0", "gene1"],
        }
    )

    custom_weight = 2.5
    custom_sign = -1.0

    graph = DIAGVI.construct_custom_guidance_graph(
        input_dict={"diss": adata_diss, "spatial": adata_spatial},
        mapping_df=mapping_df,
        weight=custom_weight,
        sign=custom_sign,
    )

    # Check that cross-modality edges have custom weight and sign
    # Self-loops have weight=1.0 and sign=1.0
    assert custom_weight in graph.edge_weight.tolist()
    assert custom_sign in graph.edge_sign.tolist()


def test_construct_custom_guidance_graph_missing_features(guidance_graph_adatas):
    """Test that missing features in mapping are handled gracefully."""
    adata_diss, adata_spatial = guidance_graph_adatas
    DIAGVI.setup_anndata(adata_diss, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nb")

    # Include a gene that doesn't exist in the data
    mapping_df = pd.DataFrame(
        {
            "diss": ["gene0", "gene1", "nonexistent_gene"],
            "spatial": ["gene0", "gene1", "gene2"],
        }
    )

    graph = DIAGVI.construct_custom_guidance_graph(
        input_dict={"diss": adata_diss, "spatial": adata_spatial},
        mapping_df=mapping_df,
    )

    # Graph should still be created, missing features are skipped
    assert graph is not None
    assert hasattr(graph, "edge_index")

    # Should only have 2 valid mapped pairs (gene0 and gene1)
    n_nodes = 20 + 15
    n_valid_mapped_pairs = 2
    expected_edges = 2 * n_valid_mapped_pairs + n_nodes
    assert graph.edge_index.shape[1] == expected_edges


def test_construct_custom_guidance_graph_feature_renaming(guidance_graph_adatas):
    """Test that features are renamed with modality suffix."""
    adata_diss, adata_spatial = guidance_graph_adatas
    DIAGVI.setup_anndata(adata_diss, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nb")

    mapping_df = pd.DataFrame(
        {
            "diss": ["gene0"],
            "spatial": ["gene0"],
        }
    )

    graph = DIAGVI.construct_custom_guidance_graph(
        input_dict={"diss": adata_diss, "spatial": adata_spatial},
        mapping_df=mapping_df,
    )

    # Check that modality indices are assigned correctly
    # First 20 nodes (diss features) should have modality_index=0
    # Next 15 nodes (spatial features) should have modality_index=1
    assert len(graph.diss_indices) == 20
    assert len(graph.spatial_indices) == 15
    assert graph.diss_indices.min() == 0
    assert graph.diss_indices.max() == 19
    assert graph.spatial_indices.min() == 20
    assert graph.spatial_indices.max() == 34


def test_construct_custom_guidance_graph_self_loops(guidance_graph_adatas):
    """Test that self-loops are created for all nodes."""
    adata_diss, adata_spatial = guidance_graph_adatas
    DIAGVI.setup_anndata(adata_diss, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nb")

    # Empty mapping - only self-loops should exist
    mapping_df = pd.DataFrame(
        {
            "diss": [],
            "spatial": [],
        }
    )

    graph = DIAGVI.construct_custom_guidance_graph(
        input_dict={"diss": adata_diss, "spatial": adata_spatial},
        mapping_df=mapping_df,
    )

    n_nodes = 20 + 15
    # With no mapped pairs, only self-loops should exist
    assert graph.edge_index.shape[1] == n_nodes

    # Check that all edges are self-loops (src == dst)
    src, dst = graph.edge_index
    assert (src == dst).all()


# =============================================================================
# Tests for DiagTrainingPlan lam_class initialization and loss_annealing
# =============================================================================


class MockModule(torch.nn.Module):
    """Mock module for testing DiagTrainingPlan initialization."""

    def __init__(self, semi_supervised):
        super().__init__()
        self.semi_supervised = semi_supervised
        # TrainingPlan expects these attributes
        self.linear = torch.nn.Linear(1, 1)

    def loss(self, *args, **kwargs):
        """Mock loss method required by TrainingPlan."""
        pass


@pytest.mark.parametrize("explicit_lam_class", [0.0, 50.0, 100.0, 200.0])
def test_lam_class_initialization_explicit(explicit_lam_class):
    """Test that explicit lam_class values are respected."""
    from scvi.external.diagvi._task import DiagTrainingPlan

    # Even with semi_supervised=True, explicit value should be used
    mock_module = MockModule(semi_supervised=True)
    plan = DiagTrainingPlan(mock_module, lam_class=explicit_lam_class)

    assert plan.lam_class == explicit_lam_class, (
        f"Expected explicit lam_class={explicit_lam_class}, got {plan.lam_class}"
    )


def test_diagvi_training_with_loss_annealing(adata_seq, adata_spatial):
    """Test training with loss_annealing=True to cover the annealing branch."""
    DIAGVI.setup_anndata(adata_seq, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nb")
    model = DIAGVI({"diss": adata_seq, "spatial": adata_spatial})

    # Train with loss_annealing enabled
    model.train(
        max_epochs=2,
        batch_size=16,
        plan_kwargs={
            "loss_annealing": True,
            "lam_sinkhorn": 1.0,
            "sinkhorn_blur": 1.0,
            "sinkhorn_reach": 1.0,
        },
    )

    assert model.is_trained_ is True


def test_anneal_param_function():
    """Test the _anneal_param helper function directly."""
    from scvi.external.diagvi._task import _anneal_param

    # Test at epoch 0 (start of training)
    result = _anneal_param(current_epoch=0, max_epochs=30, init_value=10.0, target_value=1.0)
    assert result == 10.0, f"Expected 10.0 at epoch 0, got {result}"

    # Test at epoch 5 (halfway through annealing period of 10 epochs)
    result = _anneal_param(current_epoch=5, max_epochs=30, init_value=10.0, target_value=1.0)
    assert result == 5.5, f"Expected 5.5 at epoch 5, got {result}"

    # Test at epoch 10 (end of annealing period)
    result = _anneal_param(current_epoch=10, max_epochs=30, init_value=10.0, target_value=1.0)
    assert result == 1.0, f"Expected 1.0 at epoch 10, got {result}"

    # Test at epoch 20 (after annealing period)
    result = _anneal_param(current_epoch=20, max_epochs=30, init_value=10.0, target_value=1.0)
    assert result == 1.0, f"Expected 1.0 at epoch 20, got {result}"


# =============================================================================
# Tests for different likelihoods
# =============================================================================


@pytest.fixture(scope="module")
def adata_continuous():
    """AnnData with continuous (non-negative) data for testing continuous likelihoods."""
    # Use log1p-normal data (non-negative)
    X = np.abs(np.random.randn(N_OBS_SEQ, N_VARS)) + 0.1  # Ensure positive
    obs = {"batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=N_OBS_SEQ))}
    var = pd.DataFrame(index=[f"gene{i}" for i in range(N_VARS)])
    return AnnData(X=X.astype(np.float32), obs=obs, var=var)


@pytest.fixture(scope="module")
def adata_continuous_spatial():
    """AnnData with continuous (non-negative) data for spatial modality."""
    X = np.abs(np.random.randn(N_OBS_SPATIAL, N_VARS)) + 0.1
    obs = {"batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=N_OBS_SPATIAL))}
    var = pd.DataFrame(index=[f"gene{i}" for i in range(N_VARS)])
    return AnnData(X=X.astype(np.float32), obs=obs, var=var)


@pytest.mark.parametrize("likelihood", ["nb", "zinb"])
def test_diagvi_count_likelihoods(adata_seq, adata_spatial, likelihood):
    """Test DIAGVI with different count-based likelihoods."""
    DIAGVI.setup_anndata(adata_seq, batch_key="batch", likelihood=likelihood)
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood=likelihood)
    
    model = DIAGVI({"diss": adata_seq, "spatial": adata_spatial})
    model.train(max_epochs=1, batch_size=16)
    
    assert model.is_trained_ is True
    assert model.module.modality_likelihoods["diss"] == likelihood
    assert model.module.modality_likelihoods["spatial"] == likelihood
    
    # Test latent representation
    latents = model.get_latent_representation()
    assert latents["diss"].shape[0] == adata_seq.shape[0]


@pytest.mark.parametrize("likelihood", ["log1pnormal", "ziln", "zig"])
def test_diagvi_continuous_likelihoods(adata_continuous, adata_continuous_spatial, likelihood):
    """Test DIAGVI with different continuous likelihoods."""
    DIAGVI.setup_anndata(adata_continuous, batch_key="batch", likelihood=likelihood)
    DIAGVI.setup_anndata(adata_continuous_spatial, batch_key="batch", likelihood=likelihood)
    
    model = DIAGVI({"diss": adata_continuous, "spatial": adata_continuous_spatial})
    model.train(max_epochs=1, batch_size=16)
    
    assert model.is_trained_ is True
    assert model.module.modality_likelihoods["diss"] == likelihood
    
    # Verify decoder normalization is disabled for continuous data
    assert model.module.decoder_0.normalize is False
    assert model.module.decoder_1.normalize is False


def test_diagvi_normal_likelihood(adata_continuous, adata_continuous_spatial):
    """Test DIAGVI with normal likelihood (allows negative values)."""
    DIAGVI.setup_anndata(
        adata_continuous, batch_key="batch", likelihood="normal", normalize_lib=False
    )
    DIAGVI.setup_anndata(
        adata_continuous_spatial, batch_key="batch", likelihood="normal", normalize_lib=False
    )
    
    model = DIAGVI({"diss": adata_continuous, "spatial": adata_continuous_spatial})
    model.train(max_epochs=1, batch_size=16)
    
    assert model.is_trained_ is True


def test_diagvi_nbmixture_likelihood(adata_seq, adata_spatial):
    """Test DIAGVI with nbmixture likelihood (protein-like data)."""
    DIAGVI.setup_anndata(adata_seq, batch_key="batch", likelihood="nbmixture")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nbmixture")
    
    model = DIAGVI({"diss": adata_seq, "spatial": adata_spatial})
    
    # Check dual-pathway decoder is used
    from scvi.external.diagvi._base_components import DecoderDualPathway
    assert isinstance(model.module.decoder_0, DecoderDualPathway)
    assert isinstance(model.module.decoder_1, DecoderDualPathway)
    
    model.train(max_epochs=1, batch_size=16)
    assert model.is_trained_ is True


def test_diagvi_mixed_likelihoods(adata_seq, adata_continuous_spatial):
    """Test DIAGVI with different likelihoods per modality."""
    DIAGVI.setup_anndata(adata_seq, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_continuous_spatial, batch_key="batch", likelihood="log1pnormal")
    
    model = DIAGVI({"diss": adata_seq, "spatial": adata_continuous_spatial})
    
    # Check different decoder configurations
    assert model.module.decoder_0.normalize is True  # nb uses softmax
    assert model.module.decoder_1.normalize is False  # log1pnormal does not
    
    model.train(max_epochs=1, batch_size=16)
    assert model.is_trained_ is True


# =============================================================================
# Tests for posterior_predictive_sample
# =============================================================================


def test_posterior_predictive_sample_single(trained_model):
    """Test posterior_predictive_sample with n_samples=1."""
    model, adata_seq, adata_spatial = trained_model
    
    samples = model.posterior_predictive_sample(n_samples=1)
    
    assert "diss" in samples
    assert "spatial" in samples
    assert samples["diss"].shape == (N_OBS_SEQ, N_VARS)
    assert samples["spatial"].shape == (N_OBS_SPATIAL, N_VARS)


def test_posterior_predictive_sample_multiple(trained_model):
    """Test posterior_predictive_sample with n_samples>1."""
    model, adata_seq, adata_spatial = trained_model
    n_samples = 5
    
    samples = model.posterior_predictive_sample(n_samples=n_samples)
    
    assert samples["diss"].shape == (N_OBS_SEQ, N_VARS, n_samples)
    assert samples["spatial"].shape == (N_OBS_SPATIAL, N_VARS, n_samples)


def test_posterior_predictive_sample_with_indices(trained_model):
    """Test posterior_predictive_sample with subset of indices."""
    model, adata_seq, adata_spatial = trained_model
    
    indices = {
        "diss": list(range(10)),
        "spatial": list(range(20)),
    }
    
    samples = model.posterior_predictive_sample(indices=indices, n_samples=1)
    
    assert samples["diss"].shape == (10, N_VARS)
    assert samples["spatial"].shape == (20, N_VARS)


def test_posterior_predictive_sample_invalid_modality(trained_model):
    """Test that invalid modality names raise ValueError."""
    model, adata_seq, adata_spatial = trained_model
    
    # Create fake adatas dict with invalid key
    invalid_adatas = {"invalid_modality": adata_seq}
    
    with pytest.raises(ValueError, match="Invalid modality names"):
        model.posterior_predictive_sample(adatas=invalid_adatas)


# =============================================================================
# Tests for adaptive Sinkhorn parameters
# =============================================================================


def test_diagvi_adaptive_sinkhorn(adata_seq, adata_spatial):
    """Test training with adaptive Sinkhorn parameters (blur=None, reach=None)."""
    DIAGVI.setup_anndata(adata_seq, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nb")
    model = DIAGVI({"diss": adata_seq, "spatial": adata_spatial})
    
    # Train with adaptive parameters (default: blur=None, reach=None)
    model.train(
        max_epochs=2,
        batch_size=16,
        plan_kwargs={
            "sinkhorn_blur": None,
            "sinkhorn_reach": None,
            "epsilon_from_cost": "mean",
            "epsilon_scale": 0.5,
        },
    )
    
    assert model.is_trained_ is True


@pytest.mark.parametrize("epsilon_from_cost", ["mean", "std"])
def test_diagvi_adaptive_epsilon_methods(adata_seq, adata_spatial, epsilon_from_cost):
    """Test different methods for computing adaptive epsilon."""
    DIAGVI.setup_anndata(adata_seq, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nb")
    model = DIAGVI({"diss": adata_seq, "spatial": adata_spatial})
    
    model.train(
        max_epochs=1,
        batch_size=16,
        plan_kwargs={
            "sinkhorn_blur": None,
            "epsilon_from_cost": epsilon_from_cost,
        },
    )
    
    assert model.is_trained_ is True


# =============================================================================
# Tests for deterministic inference
# =============================================================================


def test_get_imputed_values_deterministic(trained_model):
    """Test get_imputed_values with deterministic=True vs False."""
    model, adata_seq, adata_spatial = trained_model
    
    # Deterministic should give same result on repeated calls
    imputed_det1 = model.get_imputed_values(query_name="spatial", deterministic=True)
    imputed_det2 = model.get_imputed_values(query_name="spatial", deterministic=True)
    
    np.testing.assert_array_equal(imputed_det1, imputed_det2)


def test_get_latent_deterministic_mode(trained_model):
    """Test inference with deterministic mode."""
    model, adata_seq, adata_spatial = trained_model
    
    device = next(model.module.parameters()).device
    x = torch.tensor(adata_seq.X, dtype=torch.float32, device=device)
    
    # Deterministic inference should use mean instead of sampling
    outputs_det = model.module.inference(x, mode="diss", deterministic=True)
    outputs_det2 = model.module.inference(x, mode="diss", deterministic=True)
    
    # Results should be identical
    torch.testing.assert_close(outputs_det["z"], outputs_det2["z"])


# =============================================================================
# Tests for semi-supervised with both modalities
# =============================================================================


def test_semi_supervised_both_modalities(adata_seq_with_labels, adata_spatial_with_labels):
    """Test semi-supervised learning enabled for both modalities."""
    DIAGVI.setup_anndata(
        adata_seq_with_labels,
        batch_key="batch",
        labels_key="cell_type",
        likelihood="nb",
        gmm_prior=True,
    )
    DIAGVI.setup_anndata(
        adata_spatial_with_labels,
        batch_key="batch",
        labels_key="cell_type",
        likelihood="nb",
        gmm_prior=True,
    )
    
    model = DIAGVI({"diss": adata_seq_with_labels, "spatial": adata_spatial_with_labels})
    
    # Both modalities should have classifiers
    assert model.module.classifier_0 is not None
    assert model.module.classifier_1 is not None
    
    # Both should use GMM prior
    assert model.module.use_gmm_prior["diss"] is True
    assert model.module.use_gmm_prior["spatial"] is True
    
    model.train(max_epochs=1, batch_size=16)
    assert model.is_trained_ is True


# =============================================================================
# Tests for pre-provided guidance graph
# =============================================================================


def test_diagvi_with_precomputed_guidance_graph(adata_seq, adata_spatial):
    """Test DIAGVI initialization with a pre-computed guidance graph."""
    DIAGVI.setup_anndata(adata_seq, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nb")
    
    # Create a custom guidance graph
    mapping_df = pd.DataFrame({
        "diss": ["gene0", "gene1", "gene2"],
        "spatial": ["gene0", "gene1", "gene2"],
    })
    custom_graph = DIAGVI.construct_custom_guidance_graph(
        input_dict={"diss": adata_seq, "spatial": adata_spatial},
        mapping_df=mapping_df,
    )
    
    # Initialize model with pre-computed graph
    model = DIAGVI(
        {"diss": adata_seq, "spatial": adata_spatial},
        guidance_graph=custom_graph,
    )
    
    assert model.guidance_graph is custom_graph
    model.train(max_epochs=1, batch_size=16)
    assert model.is_trained_ is True


# =============================================================================
# Tests for error cases
# =============================================================================


def test_get_imputed_values_invalid_query_name(trained_model):
    """Test that invalid query_name raises ValueError."""
    model, adata_seq, adata_spatial = trained_model
    
    with pytest.raises(ValueError, match="must be one of"):
        model.get_imputed_values(query_name="invalid_modality")


def test_setup_anndata_negative_data_continuous_likelihood():
    """Test that negative data with non-normal continuous likelihood raises error."""
    # Data with negative values
    X = np.random.randn(50, 20)  # Contains negatives
    adata = AnnData(
        X=X.astype(np.float32),
        obs={"batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=50))},
    )
    
    with pytest.raises(ValueError, match="requires non-negative data"):
        DIAGVI.setup_anndata(adata, batch_key="batch", likelihood="log1pnormal")


def test_guidance_graph_consistency_check_fails():
    """Test that inconsistent guidance graph raises ValueError."""
    from scvi.external.diagvi._utils import _check_guidance_graph_consistency
    from torch_geometric.data import Data
    
    # Create adatas
    adata1 = AnnData(X=np.random.poisson(1.0, size=(50, 20)))
    adata2 = AnnData(X=np.random.poisson(1.0, size=(40, 15)))
    
    # Create graph with wrong number of nodes
    wrong_graph = Data(
        x=torch.eye(10),  # Wrong: should be 20+15=35
        edge_index=torch.tensor([[0, 1], [1, 0]]).T,
        edge_weight=torch.ones(2),
        edge_sign=torch.ones(2),
    )
    
    with pytest.raises(ValueError, match="node count"):
        _check_guidance_graph_consistency(wrong_graph, {"mod1": adata1, "mod2": adata2})


def test_guidance_graph_missing_self_loops():
    """Test that graph without self-loops raises ValueError."""
    from scvi.external.diagvi._utils import _check_guidance_graph_consistency
    from torch_geometric.data import Data
    
    adata1 = AnnData(X=np.random.poisson(1.0, size=(50, 5)))
    adata2 = AnnData(X=np.random.poisson(1.0, size=(40, 5)))
    
    # Graph with no self-loops
    no_selfloop_graph = Data(
        x=torch.eye(10),
        edge_index=torch.tensor([[0, 1], [1, 0]]).T,  # No self-loops
        edge_weight=torch.ones(2),
        edge_sign=torch.ones(2),
    )
    
    with pytest.raises(ValueError, match="self-loops"):
        _check_guidance_graph_consistency(no_selfloop_graph, {"mod1": adata1, "mod2": adata2})


# =============================================================================
# Tests for base components
# =============================================================================


def test_decoder_single_pathway_outputs():
    """Test DecoderSinglePathway output shapes."""
    from scvi.external.diagvi._base_components import DecoderSinglePathway
    
    n_output = 50
    n_batches = 2
    n_latent = 10
    batch_size = 32
    
    decoder = DecoderSinglePathway(n_output=n_output, n_batches=n_batches, normalize=True)
    
    u = torch.randn(batch_size, n_latent)
    l = torch.randn(batch_size, 1)
    batch_index = torch.randint(0, n_batches, (batch_size,))
    v = torch.randn(n_output, n_latent)
    
    px_scale, px_r, px_rate, px_dropout = decoder(u, l, batch_index, v)
    
    assert px_scale.shape == (batch_size, n_output)
    assert px_r.shape == (batch_size, n_output)
    assert px_rate.shape == (batch_size, n_output)
    assert px_dropout.shape == (n_output,)
    
    # Check softmax normalization (sums to 1)
    torch.testing.assert_close(px_scale.sum(dim=1), torch.ones(batch_size), atol=1e-5, rtol=1e-5)


def test_decoder_dual_pathway_outputs():
    """Test DecoderDualPathway output shapes."""
    from scvi.external.diagvi._base_components import DecoderDualPathway
    
    n_output = 50
    n_batches = 2
    n_latent = 10
    batch_size = 32
    
    decoder = DecoderDualPathway(n_output=n_output, n_batches=n_batches, normalize=True)
    
    u = torch.randn(batch_size, n_latent)
    l = torch.randn(batch_size, 1)
    batch_index = torch.randint(0, n_batches, (batch_size,))
    v = torch.randn(n_output, n_latent)
    
    scales, px_r, rates, mixture_logits = decoder(u, l, batch_index, v)
    
    # Dual pathway returns tuples
    assert isinstance(scales, tuple) and len(scales) == 2
    assert isinstance(rates, tuple) and len(rates) == 2
    
    assert scales[0].shape == (batch_size, n_output)
    assert scales[1].shape == (batch_size, n_output)
    assert rates[0].shape == (batch_size, n_output)
    assert rates[1].shape == (batch_size, n_output)
    assert px_r.shape == (batch_size, n_output)
    assert mixture_logits.shape == (batch_size, n_output)


def test_graph_encoder_outputs():
    """Test GraphEncoder output shapes."""
    from scvi.external.diagvi._base_components import GraphEncoder
    
    n_nodes = 100
    n_latent = 50
    
    encoder = GraphEncoder(vnum=n_nodes, out_features=n_latent)
    
    # Create simple edge index with self-loops
    edge_index = torch.stack([torch.arange(n_nodes), torch.arange(n_nodes)])
    
    z, mu, logvar = encoder(edge_index)
    
    assert z.shape == (n_nodes, n_latent)
    assert mu.shape == (n_nodes, n_latent)
    assert logvar.shape == (n_nodes, n_latent)


def test_decoder_batch_index_out_of_bounds():
    """Test that out-of-bounds batch index raises IndexError."""
    from scvi.external.diagvi._base_components import DecoderSinglePathway
    
    decoder = DecoderSinglePathway(n_output=50, n_batches=2, normalize=True)
    
    u = torch.randn(10, 10)
    l = torch.randn(10, 1)
    v = torch.randn(50, 10)
    invalid_batch = torch.tensor([0, 1, 2, 3, 0, 1, 2, 3, 0, 1])  # 2 and 3 are out of bounds
    
    with pytest.raises(IndexError, match="out of bounds"):
        decoder(u, l, invalid_batch, v)


# =============================================================================
# Tests for utility functions
# =============================================================================


def test_construct_guidance_graph_auto():
    """Test automatic guidance graph construction with shared features."""
    from scvi.external.diagvi._utils import _construct_guidance_graph
    
    # Create adatas with overlapping features
    shared_genes = [f"shared_gene{i}" for i in range(10)]
    unique_genes1 = [f"unique1_gene{i}" for i in range(5)]
    unique_genes2 = [f"unique2_gene{i}" for i in range(5)]
    
    adata1 = AnnData(
        X=np.random.poisson(1.0, size=(50, 15)),
        var=pd.DataFrame(index=shared_genes + unique_genes1),
    )
    adata2 = AnnData(
        X=np.random.poisson(1.0, size=(40, 15)),
        var=pd.DataFrame(index=shared_genes + unique_genes2),
    )
    
    graph = _construct_guidance_graph({"mod1": adata1, "mod2": adata2}, mapping_df=None)
    
    # Should have 15+15=30 nodes
    assert graph.num_nodes == 30
    
    # Edges: 10 shared * 2 (bidirectional) + 30 self-loops = 50
    assert graph.edge_index.shape[1] == 50


def test_construct_guidance_graph_no_overlap_error():
    """Test that no overlapping features raises ValueError."""
    from scvi.external.diagvi._utils import _construct_guidance_graph
    
    adata1 = AnnData(
        X=np.random.poisson(1.0, size=(50, 5)),
        var=pd.DataFrame(index=[f"gene1_{i}" for i in range(5)]),
    )
    adata2 = AnnData(
        X=np.random.poisson(1.0, size=(40, 5)),
        var=pd.DataFrame(index=[f"gene2_{i}" for i in range(5)]),  # No overlap
    )
    
    with pytest.raises(ValueError, match="No overlapping features"):
        _construct_guidance_graph({"mod1": adata1, "mod2": adata2}, mapping_df=None)


def test_kl_divergence_graph():
    """Test KL divergence computation for graph latent variables."""
    from scvi.external.diagvi._utils import kl_divergence_graph
    
    n_nodes = 100
    n_latent = 50
    
    # Standard normal: mu=0, logvar=0 -> KL should be 0
    mu_zero = torch.zeros(n_nodes, n_latent)
    logvar_zero = torch.zeros(n_nodes, n_latent)
    
    kl_zero = kl_divergence_graph(mu_zero, logvar_zero)
    assert kl_zero.item() < 1e-5, f"KL for standard normal should be ~0, got {kl_zero.item()}"
    
    # Non-zero mu and logvar should give positive KL
    mu_nonzero = torch.randn(n_nodes, n_latent)
    logvar_nonzero = torch.randn(n_nodes, n_latent)
    
    kl_nonzero = kl_divergence_graph(mu_nonzero, logvar_nonzero)
    assert kl_nonzero.item() > 0, "KL should be positive for non-standard distribution"
