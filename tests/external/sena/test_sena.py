import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.external import SENADVAE


@pytest.fixture
def mock_sena_adata():
    """Create mock AnnData with perturbation annotations for SENADVAE testing."""
    adata = synthetic_iid(
        n_genes=100,
    )

    # Add perturbation annotations
    # Create a mix of control and perturbed cells
    n_cells = adata.n_obs
    perturbations = [""] * (n_cells // 2)  # Half control cells

    # Add some single perturbations
    gene_names = ["GENE1", "GENE2", "GENE3"]
    for i in range(n_cells // 2, n_cells):
        perturbations.append(gene_names[i % len(gene_names)])

    adata.obs["perturbation"] = perturbations

    return adata


@pytest.fixture
def mock_go_files(tmp_path):
    """Create mock GO annotation files for SENADVAE testing."""
    # Create mock GO pathway file
    go_file_path = tmp_path / "go_pathways.tsv"
    go_pathways = ["GO:0001", "GO:0002", "GO:0003"]
    with open(go_file_path, "w") as f:
        f.write("PathwayID\n")
        for go in go_pathways:
            f.write(f"{go}\n")

    # Create mock gene-to-GO mapping file
    go_gene_map_path = tmp_path / "gene_go_map.tsv"
    with open(go_gene_map_path, "w") as f:
        f.write("GO_id\tensembl_id\n")
        # Map some genes to GO terms (using var_names from synthetic data)
        for i, go in enumerate(go_pathways):
            for j in range(10):
                gene_idx = (i * 10 + j) % 100
                f.write(f"{go}\tgene_{gene_idx}\n")

    # Create mock gene symbol to ensemble mapping file
    gene_symb_ensemble_path = tmp_path / "gene_symb_ensemble.tsv"
    with open(gene_symb_ensemble_path, "w") as f:
        f.write("external_gene_name\tensembl_gene_id\n")
        f.write("GENE1\tgene_0\n")
        f.write("GENE2\tgene_10\n")
        f.write("GENE3\tgene_20\n")

    return {
        "go_file_path": str(go_file_path),
        "go_gene_map_path": str(go_gene_map_path),
        "gene_symb_ensemble_path": str(gene_symb_ensemble_path),
    }


def test_senadvae_setup_anndata(mock_sena_adata):
    """Test SENADVAE.setup_anndata registers data correctly."""
    SENADVAE.setup_anndata(mock_sena_adata, perturbation_key="perturbation")

    assert "_sena_perturbation_key" in mock_sena_adata.uns
    assert mock_sena_adata.uns["_sena_perturbation_key"] == "perturbation"
    assert "intervention_genes" in mock_sena_adata.uns
    assert "n_intervention_genes" in mock_sena_adata.uns


def test_senadvae_analyze_perturbations(mock_sena_adata):
    """Test SENADVAE.analyze_perturbations provides correct statistics."""
    stats = SENADVAE.analyze_perturbations(mock_sena_adata, "perturbation")

    assert "n_total_cells" in stats
    assert "n_controls" in stats
    assert "n_unique_genes" in stats
    assert stats["n_total_cells"] == mock_sena_adata.n_obs
    assert stats["n_controls"] > 0


def test_senadvae_init(mock_sena_adata, mock_go_files):
    """Test SENADVAE model initialization."""
    SENADVAE.setup_anndata(mock_sena_adata, perturbation_key="perturbation")

    model = SENADVAE(
        mock_sena_adata,
        go_file_path=mock_go_files["go_file_path"],
        go_gene_map_path=mock_go_files["go_gene_map_path"],
        gene_symb_ensemble_path=mock_go_files["gene_symb_ensemble_path"],
        n_hidden_encoder=32,
        n_hidden_decoder=32,
        n_hidden_interv=32,
        n_go_thresh=2,
        seed=42,
    )

    assert model is not None
    assert model.module is not None


def test_senadvae_train(mock_sena_adata, mock_go_files):
    """Test SENADVAE model training."""
    SENADVAE.setup_anndata(mock_sena_adata, perturbation_key="perturbation")

    model = SENADVAE(
        mock_sena_adata,
        go_file_path=mock_go_files["go_file_path"],
        go_gene_map_path=mock_go_files["go_gene_map_path"],
        gene_symb_ensemble_path=mock_go_files["gene_symb_ensemble_path"],
        n_hidden_encoder=32,
        n_hidden_decoder=32,
        n_hidden_interv=32,
        n_go_thresh=2,
        seed=42,
    )

    model.train(max_epochs=2, batch_size=32, check_val_every_n_epoch=1)

    assert model.is_trained_


def test_senadvae_get_latent_representation(mock_sena_adata, mock_go_files):
    """Test SENADVAE latent representation extraction."""
    SENADVAE.setup_anndata(mock_sena_adata, perturbation_key="perturbation")

    model = SENADVAE(
        mock_sena_adata,
        go_file_path=mock_go_files["go_file_path"],
        go_gene_map_path=mock_go_files["go_gene_map_path"],
        gene_symb_ensemble_path=mock_go_files["gene_symb_ensemble_path"],
        n_hidden_encoder=32,
        n_hidden_decoder=32,
        n_hidden_interv=32,
        n_go_thresh=2,
        seed=42,
    )

    model.train(max_epochs=1, batch_size=32)

    latent = model.get_latent_representation()

    assert latent is not None
    assert isinstance(latent, np.ndarray)
    assert latent.shape[0] == mock_sena_adata.n_obs // 2


def test_senadvae_get_causal_graph(mock_sena_adata, mock_go_files):
    """Test SENADVAE causal graph extraction."""
    SENADVAE.setup_anndata(mock_sena_adata, perturbation_key="perturbation")

    model = SENADVAE(
        mock_sena_adata,
        go_file_path=mock_go_files["go_file_path"],
        go_gene_map_path=mock_go_files["go_gene_map_path"],
        gene_symb_ensemble_path=mock_go_files["gene_symb_ensemble_path"],
        n_hidden_encoder=32,
        n_hidden_decoder=32,
        n_hidden_interv=32,
        n_go_thresh=2,
        seed=42,
    )

    model.train(max_epochs=1, batch_size=32)

    causal_graph = model.get_causal_graph()

    assert causal_graph is not None
    assert isinstance(causal_graph, np.ndarray)
