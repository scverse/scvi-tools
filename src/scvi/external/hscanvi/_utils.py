import pronto
import networkx as nx
import os
import requests
import numpy as np
from tqdm import tqdm


class CellOntologyNavigator:
    def __init__(self, cts_of_interest: list[str]):
        self.ontology_url = "http://purl.obolibrary.org/obo/cl.owl"
        self.ontology_path = self.download_ontology_if_necessary()
        self.full_ontology = self._load_cell_ontology_as_dag()

        self.nodes_of_interest = self._obtain_valid_nodes(cts_of_interest)
        self.multilabel_matrix = self.get_multilabel_matrix()

        self.node_to_index = {
            node: idx for idx, node in enumerate(self.nodes_of_interest)
        }

        name_to_node = {}
        for node in self.full_ontology.nodes():
            node_label = self.full_ontology.nodes[node]["name"]  # assumes no duplicates
            name_to_node[node_label] = node
        self.name_to_node = name_to_node

    def download_ontology_if_necessary(self) -> str:
        local_filename = "cl.owl"
        if not os.path.exists(local_filename):
            print(f"Downloading ontology from {self.ontology_url}...")
            try:
                response = requests.get(self.ontology_url, stream=True)
                response.raise_for_status()
                with open(local_filename, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                print("Download complete.")
            except requests.exceptions.RequestException as e:
                print(f"Error downloading ontology: {e}")
                return None
        else:
            print("Ontology file already exists locally.")
        return local_filename

    def _load_cell_ontology_as_dag(self) -> nx.DiGraph:
        try:
            ontology = pronto.Ontology(self.ontology_path)

            # build graph containing only CL terms
            graph = nx.DiGraph()
            for term in ontology.terms():
                stop = (not term.id.startswith("CL:")) or term.name is None
                if stop:
                    continue
                graph.add_node(term.id, name=term.name, definition=term.definition)
                for superclass in term.superclasses():
                    if not superclass.id.startswith("CL:"):
                        continue
                    if superclass.id in ontology:
                        if term.id == superclass.id:
                            continue
                        graph.add_edge(term.id, superclass.id, type="is_a")

            print(
                f"Successfully loaded ontology with {len(graph.nodes)} nodes and {len(graph.edges)} edges."
            )
            if not nx.is_directed_acyclic_graph(graph):
                raise ValueError("The ontology graph is not a DAG.")
            else:
                print("The ontology graph is a DAG; resuming...")
            return graph

        except Exception as e:
            print(f"An error occurred: {e}")
            return None

    def _obtain_valid_nodes(self, cts_of_interest: list[str]) -> list[str]:
        valid_nodes = []
        ontology = self.full_ontology.copy()
        for node in ontology.nodes():
            node_name = ontology.nodes[node]["name"]
            if node_name in cts_of_interest:
                valid_nodes.append(node)

        nodes_to_include = set(valid_nodes)
        for node in valid_nodes:
            nodes_to_include.update(nx.ancestors(ontology, node))

        nodes_to_include = list(nodes_to_include)
        return nodes_to_include

    def get_multilabel_matrix(self) -> np.ndarray:
        mat = np.zeros((len(self.nodes_of_interest), len(self.nodes_of_interest)))
        for idx_a, node_a in enumerate(tqdm(self.nodes_of_interest)):
            for idx_b, node_b in enumerate(self.nodes_of_interest):
                mat[idx_a, idx_b] = float(
                    nx.has_path(self.full_ontology, node_a, node_b)
                )
        return mat

    def populate_adata(
        self,
        adata,
        cell_type_key: str = "cell_type",
        added_node_key: str = "hscanvi_co_node",
        added_node_idx_key: str = "hscanvi_co_node_idx",
        added_multilabel_key: str = "hscanvi_co_multilabel_onehot",
    ):
        if cell_type_key not in adata.obs:
            raise ValueError(f"{cell_type_key} not present in `adata.obs`.")
        adata.obs.loc[:, added_node_key] = adata.obs[cell_type_key].map(
            self.name_to_node
        )
        adata.obs.loc[:, added_node_idx_key] = adata.obs[added_node_key].map(
            self.node_to_index
        )

        multilabel_onehot = self.multilabel_matrix[
            adata.obs[added_node_idx_key].values, :
        ]
        adata.obsm[added_multilabel_key] = multilabel_onehot


def process_adata_ontology(
    adata,
    cell_type_key: str = "cell_type",
    added_node_key: str = "hscanvi_co_node",
    added_node_idx_key: str = "hscanvi_co_node_idx",
    added_multilabel_key: str = "hscanvi_co_multilabel_onehot",
):
    cell_ontology = CellOntologyNavigator(adata.obs[cell_type_key].unique())
    cell_ontology.populate_adata(
        adata, cell_type_key, added_node_key, added_node_idx_key, added_multilabel_key
    )
