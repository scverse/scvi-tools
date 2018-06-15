from .brain_large import BrainLargeDataset
from .cortex import CortexDataset
from .dataset import GeneExpressionDataset
from .synthetic import SyntheticDataset
from .retina import RetinaDataset
from .cite_seq import CiteSeqDataset
from .brain_small import BrainSmallDataset
from .hemato import HematoDataset
from .pbmc import PbmcDataset
from .loom import LoomDataset

__all__ = ['SyntheticDataset',
           'CortexDataset',
           'BrainLargeDataset',
           'RetinaDataset',
           'GeneExpressionDataset',
           'CiteSeqDataset',
           'BrainSmallDataset',
           'HematoDataset',
           'PbmcDataset',
           'LoomDataset']


def load_datasets(dataset_name, save_path='data/'):
    if dataset_name == 'synthetic':
        gene_dataset = SyntheticDataset()
    elif dataset_name == 'cortex':
        gene_dataset = CortexDataset()
    elif dataset_name == 'brain_large':
        gene_dataset = BrainLargeDataset(subsample_size=128, save_path=save_path)
    elif dataset_name == 'retina':
        gene_dataset = RetinaDataset(save_path=save_path)
    elif dataset_name == 'cite_seq_cbmc' or dataset_name == 'cite_seq_pbmc':
        gene_dataset = CiteSeqDataset(name=dataset_name.split('_')[-1], save_path=save_path)
    elif dataset_name == 'brain_small':
        gene_dataset = BrainSmallDataset(save_path=save_path)
    elif dataset_name == 'hemato':
        gene_dataset = HematoDataset(save_path=save_path)
    elif dataset_name == 'pbmc':
        gene_dataset = PbmcDataset(save_path=save_path)
    elif dataset_name[-5:] == ".loom":
        gene_dataset = LoomDataset(filename=dataset_name, save_path=save_path)
    else:
        raise "No such dataset available"
    return gene_dataset
