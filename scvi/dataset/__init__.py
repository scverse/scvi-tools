from .brain_large import BrainLargeDataset
from .cortex import CortexDataset
from .dataset import GeneExpressionDataset
from .synthetic import SyntheticDataset

__all__ = ['SyntheticDataset',
           'CortexDataset',
           'BrainLargeDataset',
           'GeneExpressionDataset']


def load_datasets(dataset_name, test=False):
    if dataset_name == 'synthetic':
        gene_dataset_train, gene_dataset_test = SyntheticDataset(), SyntheticDataset()
    elif dataset_name == 'cortex':
        gene_dataset_train, gene_dataset_test = CortexDataset(type="train"), CortexDataset(type="test")
    elif dataset_name == 'brain_large':
        gene_dataset = BrainLargeDataset(subsample_size=128, test=test)
        gene_dataset_train, gene_dataset_test = gene_dataset, gene_dataset
    else:
        raise "No such dataset available"
    return gene_dataset_train, gene_dataset_test
