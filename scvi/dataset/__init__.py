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
        gene_dataset_train, gene_dataset_test = SyntheticDataset.get_dataset(), SyntheticDataset.get_dataset()
    elif dataset_name == 'cortex':
        gene_dataset_train = CortexDataset.get_dataset(type="train")
        gene_dataset_test = CortexDataset.get_dataset(type="test")
    elif dataset_name == 'brain_large':
        gene_dataset = BrainLargeDataset.get_dataset(subsample_size=1)
        gene_dataset_train, gene_dataset_test = gene_dataset, gene_dataset
    else:
        raise "No such dataset available"
    return gene_dataset_train, gene_dataset_test
