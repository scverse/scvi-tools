import numpy as np

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
    elif dataset_name.endswith('.npy'):
        data = np.load(dataset_name)
        train_data, test_data = GeneExpressionDataset.train_test_split(data)
        gene_dataset_train, gene_dataset_test = GeneExpressionDataset([train_data]), GeneExpressionDataset([test_data])
    else:
        raise "No such dataset available"
    return gene_dataset_train, gene_dataset_test
