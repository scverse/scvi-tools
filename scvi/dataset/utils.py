import numpy as np
from torch.utils.data import DataLoader
from torch.utils.data.sampler import SubsetRandomSampler

from scvi.dataset import GeneExpressionDataset


def concat_datasets(*gene_datasets):
    """
    Combines multiple unlabelled gene_datasets based on the intersection of gene names intersection.
    Datasets should all have gene_dataset.n_labels=0.
    Batch indices are generated in the same order as datasets are given.
    :param gene_datasets: a sequence of gene_datasets object
    :return: a GeneExpressionDataset instance of the concatenated datasets
    """
    assert all([hasattr(gene_dataset, 'gene_names') for gene_dataset in gene_datasets])

    def filter_genes(gene_dataset, gene_names_ref):
        """
        :return: gene_dataset.X filtered by the corresponding genes ( / columns / features)
        """
        gene_names = list(gene_dataset.gene_names)
        idx_genes = np.array([gene_names.index(gene_name) for gene_name in gene_names_ref], dtype=np.byte)
        print(len(idx_genes))
        return gene_dataset.X[:, idx_genes]

    gene_names = set.intersection(*[set(gene_dataset.gene_names) for gene_dataset in gene_datasets])
    GeneExpressionDataset(
        *GeneExpressionDataset.get_attributes_from_list(
            [filter_genes(gene_dataset, gene_names) for gene_dataset in gene_datasets]
        ), gene_names=np.array(list(gene_names)))


def get_indices(labels, n_labelled_samples_per_class_array):
    labels = np.array(labels).ravel()
    np.random.seed(0)
    permutation_idx = np.random.permutation(len(labels))
    labels = labels[permutation_idx]
    indices = []
    current_nbrs = np.zeros(len(n_labelled_samples_per_class_array))
    for idx, (label) in enumerate(labels):
        label = int(label)
        if current_nbrs[label] < n_labelled_samples_per_class_array[label]:
            indices.insert(0, idx)
            current_nbrs[label] += 1
        else:
            indices.append(idx)
    indices = np.array(indices)
    total_labelled = sum(n_labelled_samples_per_class_array)
    return permutation_idx[indices[:total_labelled]], permutation_idx[indices[total_labelled:]]


def get_data_loaders(gene_dataset, n_labelled_samples_per_class, **kwargs_data_loader):
    """
    For use in train_semi_supervised
    Get labelled/unlabelled/all data_loaders from gene_dataset
    :param gene_dataset: the dataset to consider
    :param n_labelled_samples_per_class: number of labelled samples per class
    :param kwargs_data_loader: any additional keyword arguments to pass to the dataloaders
    :return: data_loader_labelled/data_loader_unlabelled/data_loader_all
    """
    indices_labelled, indices_unlabelled = (
        get_indices(gene_dataset.labels, [n_labelled_samples_per_class] * gene_dataset.n_labels)
    )
    data_loader_all = DataLoader(gene_dataset, collate_fn=gene_dataset.collate_fn, shuffle=True, **kwargs_data_loader)
    data_loader_labelled = DataLoader(gene_dataset, collate_fn=gene_dataset.collate_fn,
                                      sampler=SubsetRandomSampler(indices_labelled), **kwargs_data_loader)
    data_loader_unlabelled = DataLoader(gene_dataset, collate_fn=gene_dataset.collate_fn,
                                        sampler=SubsetRandomSampler(indices_unlabelled), **kwargs_data_loader)
    return data_loader_all, data_loader_labelled, data_loader_unlabelled


def get_raw_data(*data_loaders):
    """
    From a list of data_loaders, return the numpy array suitable for classification with sklearn
    :param data_loaders: a sequence of data_loaders arguments
    :return: a sequence of (data, labels) for each data_loader
    """

    def get_data_labels(data_loader):
        if hasattr(data_loader, 'sampler'):
            if hasattr(data_loader.sampler, 'indices'):
                data = data_loader.dataset.X[data_loader.sampler.indices]
                labels = data_loader.dataset.labels[data_loader.sampler.indices]
        else:
            data = data_loader.dataset.X
            labels = data_loader.dataset.labels
        return data, labels.ravel()

    return [get_data_labels(data_loader) for data_loader in data_loaders]
