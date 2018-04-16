import numpy as np

from .brain_large import BrainLargeDataset
from .cortex import CortexDataset
from .dataset import GeneExpressionDataset
from .synthetic import SyntheticDataset

__all__ = ['SyntheticDataset',
           'CortexDataset',
           'BrainLargeDataset',
           'GeneExpressionDataset']


def load_datasets(dataset_name):
    if dataset_name == 'synthetic':
        gene_dataset_train, gene_dataset_test = SyntheticDataset(), SyntheticDataset()
    elif dataset_name == 'cortex':
        gene_dataset_train, gene_dataset_test = CortexDataset(type="train"), CortexDataset(type="test")
    elif dataset_name == 'brain_large':
        gene_dataset = BrainLargeDataset()
        gene_dataset_train, gene_dataset_test = gene_dataset, gene_dataset  # Return same object for now
    elif dataset_name.endswith('.npy'):
        data = np.load(dataset_name)
        train_data, test_data = GeneExpressionDataset.train_test_split(data)
        gene_dataset_train, gene_dataset_test = GeneExpressionDataset([train_data]), GeneExpressionDataset([test_data])
    else:
        raise "No such dataset available"
    return gene_dataset_train, gene_dataset_test
