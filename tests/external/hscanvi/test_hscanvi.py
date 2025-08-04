import scanpy as sc

from scvi.external.hscanvi._model import HSCANVI
from scvi.external.hscanvi._utils import CellOntologyNavigator, process_adata_ontology

path_to_adata = (
    "https://datasets.cellxgene.cziscience.com/130c04bd-8384-4efc-a529-1ed809762adc.h5ad"
)
filename = "linarsson_immune.h5ad"


def test_cell_ontology_navigator():
    adata = sc.read(filename, backup_url=path_to_adata)
    cts_of_interest = adata.obs["cell_type"].unique()

    cell_ontology = CellOntologyNavigator(cts_of_interest)
    cell_ontology.populate_adata(adata)


def test_process_adata_ontology():
    adata = sc.read(filename, backup_url=path_to_adata)
    process_adata_ontology(
        adata,
        cell_type_key="cell_type",
        added_node_key="hscanvi_co_node",
        added_node_idx_key="hscanvi_co_node_idx",
        added_multilabel_key="hscanvi_co_multilabel_onehot",
    )
    assert "hscanvi_co_node" in adata.obs
    assert "hscanvi_co_node_idx" in adata.obs
    assert "hscanvi_co_multilabel_onehot" in adata.obsm


def test_hscanvi():
    adata = sc.read(filename, backup_url=path_to_adata)
    HSCANVI.setup_anndata(adata, labels_key="cell_type", unlabeled_category="Unknown")

    model = HSCANVI(adata)
    model.train(max_epochs=2)
