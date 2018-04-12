from .cortex import CortexDataset
from .dataset import GeneExpressionDataset
from .synthetic import SyntheticDataset

__all__ = ['SyntheticDataset',
           'CortexDataset',
           'GeneExpressionDataset']


def load_datasets(dataset_name):
    if dataset_name == 'synthetic':
        gene_dataset_train, gene_dataset_test = SyntheticDataset(), SyntheticDataset()
    elif dataset_name == 'cortex':
        gene_dataset_train, gene_dataset_test = CortexDataset(type="train"), CortexDataset(type="test")
    else:
        raise "No such dataset available"
    return gene_dataset_train, gene_dataset_test
