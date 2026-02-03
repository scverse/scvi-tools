import numpy as np
import pandas as pd
import pytest
import torch
from anndata import AnnData
from mudata import MuData

from scvi.external.diagvi._model import DIAGVI

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
    loss = model.module.loss(tensors, inference_outputs, generative_outputs, mode="diss")
    assert loss.loss is not None
    assert loss.reconstruction_loss is not None
    assert loss.kl_local is not None


def test_diagvi_get_imputed_values(adata_seq, adata_spatial):
    model = make_model(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)
    # Impute spatial from diss
    imputed = model.get_imputed_values(source_name="spatial")
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape[0] == adata_spatial.shape[0]
    assert imputed.shape[1] == adata_seq.shape[1]
    # Impute diss from spatial
    imputed2 = model.get_imputed_values(source_name="diss")
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


@pytest.mark.parametrize("target_batch", [0, 1, "batch1", "batch2"])
def test_get_imputed_values_target_batch_scalar(trained_model, target_batch):
    """Test get_imputed_values with scalar target_batch (int or string)."""
    model, adata_seq, adata_spatial = trained_model

    # Impute from spatial to diss
    imputed = model.get_imputed_values(source_name="spatial", target_batch=target_batch)
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape == (N_OBS_SPATIAL, N_VARS)

    # Impute from diss to spatial
    imputed2 = model.get_imputed_values(source_name="diss", target_batch=target_batch)
    assert isinstance(imputed2, np.ndarray)
    assert imputed2.shape == (N_OBS_SEQ, N_VARS)


@pytest.mark.parametrize(
    "batch_array_type",
    ["int_array", "string_array"],
)
def test_get_imputed_values_target_batch_array(trained_model, batch_array_type):
    """Test get_imputed_values with array-like target_batch."""
    model, adata_seq, adata_spatial = trained_model

    if batch_array_type == "int_array":
        # Integer array matching spatial obs count
        target_batch = np.random.choice([0, 1], size=N_OBS_SPATIAL)
    else:
        # String array matching spatial obs count
        target_batch = np.random.choice(["batch1", "batch2"], size=N_OBS_SPATIAL)

    imputed = model.get_imputed_values(source_name="spatial", target_batch=target_batch)
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape == (N_OBS_SPATIAL, N_VARS)


@pytest.mark.parametrize("target_libsize", [100.0, 1000.0, 10000.0])
def test_get_imputed_values_target_libsize_scalar(trained_model, target_libsize):
    """Test get_imputed_values with scalar target_libsize."""
    model, adata_seq, adata_spatial = trained_model

    imputed = model.get_imputed_values(source_name="spatial", target_libsize=target_libsize)
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape == (N_OBS_SPATIAL, N_VARS)


def test_get_imputed_values_target_libsize_array(trained_model):
    """Test get_imputed_values with array target_libsize."""
    model, adata_seq, adata_spatial = trained_model

    # Array of library sizes matching obs count
    target_libsize = np.random.uniform(500, 2000, size=N_OBS_SPATIAL)

    imputed = model.get_imputed_values(source_name="spatial", target_libsize=target_libsize)
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape == (N_OBS_SPATIAL, N_VARS)


def test_get_imputed_values_combined_batch_and_libsize(trained_model):
    """Test get_imputed_values with both target_batch and target_libsize specified."""
    model, adata_seq, adata_spatial = trained_model

    # Use scalar values
    imputed = model.get_imputed_values(
        source_name="spatial", target_batch=1, target_libsize=1000.0
    )
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape == (N_OBS_SPATIAL, N_VARS)

    # Use array values
    target_batch_arr = np.random.choice([0, 1], size=N_OBS_SPATIAL)
    target_libsize_arr = np.random.uniform(500, 2000, size=N_OBS_SPATIAL)

    imputed2 = model.get_imputed_values(
        source_name="spatial",
        target_batch=target_batch_arr,
        target_libsize=target_libsize_arr,
    )
    assert isinstance(imputed2, np.ndarray)
    assert imputed2.shape == (N_OBS_SPATIAL, N_VARS)

    # Mixed: scalar batch, array libsize
    imputed3 = model.get_imputed_values(
        source_name="spatial",
        target_batch="batch1",
        target_libsize=target_libsize_arr,
    )
    assert isinstance(imputed3, np.ndarray)
    assert imputed3.shape == (N_OBS_SPATIAL, N_VARS)


@pytest.mark.parametrize(
    ("error_case", "error_kwargs"),
    [
        ("wrong_size_batch", {"target_batch": np.array([0, 1, 0])}),  # Wrong size
        ("wrong_size_libsize", {"target_libsize": np.array([100.0, 200.0, 300.0])}),  # Wrong size
        ("2d_libsize", {"target_libsize": np.array([[100.0], [200.0], [300.0]])}),  # 2D array
    ],
)
def test_get_imputed_values_errors(trained_model, error_case, error_kwargs):
    """Test that get_imputed_values raises ValueError for invalid inputs."""
    model, adata_seq, adata_spatial = trained_model

    with pytest.raises(ValueError):
        model.get_imputed_values(source_name="spatial", **error_kwargs)


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
    imputed = model.get_imputed_values(source_name="spatial")
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape[0] == N_OBS_SPATIAL

    # Impute from rna to spatial
    imputed2 = model.get_imputed_values(source_name="rna")
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
    DIAGVI.setup_anndata(
        adata_seq_with_labels,
        batch_key="batch",
        labels_key="cell_type",
        likelihood="nb",
        semi_supervised=True,
        gmm_prior=True,
    )
    DIAGVI.setup_anndata(
        adata_spatial, batch_key="batch", likelihood="nb", semi_supervised=False, gmm_prior=False
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
