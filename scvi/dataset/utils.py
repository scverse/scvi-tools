from itertools import cycle

import numpy as np
from sklearn.model_selection._split import _validate_shuffle_split
from torch.utils.data import DataLoader
from torch.utils.data.sampler import SubsetRandomSampler, SequentialSampler, RandomSampler


def filter_genes(gene_dataset, gene_names_ref, on='gene_names'):
    """
    :return: gene_dataset.X filtered by the corresponding genes ( / columns / features), idx_genes
    """
    gene_names = list(getattr(gene_dataset, on))
    subset_genes = np.array([gene_names.index(gene_name) for gene_name in gene_names_ref], dtype=np.int64)
    return gene_dataset.X[:, subset_genes], subset_genes


def arrange_categories(original_categories, mapping_from=None, mapping_to=None):
    unique_categories = np.unique(original_categories)
    n_categories = len(unique_categories)
    if mapping_to is None:
        mapping_to = range(n_categories)
    if mapping_from is None:
        mapping_from = unique_categories
    assert n_categories == len(mapping_from)
    assert n_categories == len(mapping_to)

    new_categories = np.copy(original_categories)
    for idx_from, idx_to in zip(mapping_from, mapping_to):
        new_categories[original_categories == idx_from] = idx_to
    return new_categories, n_categories


def get_indices(labels, n_labelled_samples_per_class_array, seed=0):
    labels = np.array(labels).ravel()
    np.random.seed(seed=seed)
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


def zip_first_cycle(*args):
    return zip(args[0], *[cycle(arg) for arg in args[1:]])


class DataLoaders:
    to_monitor = []
    data_loaders_loop = []

    def __init__(self, gene_dataset, **data_loaders_kwargs):
        """
        :param gene_dataset: a GeneExpressionDataset instance
        :param data_loaders_kwargs: any additional keyword arguments to pass to the data_loaders at .init
        """
        self.gene_dataset = gene_dataset
        self.data_loaders_kwargs = {
            "batch_size": 128,
            "pin_memory": True,
            "collate_fn": self.gene_dataset.collate_fn
        }
        self.data_loaders_kwargs.update(data_loaders_kwargs)
        self.data_loaders_dict = {'sequential': self()}
        self.infinite_random_sampler = cycle(self(shuffle=True))  # convenient for debugging - see .sample()

    def sample(self):
        return next(self.infinite_random_sampler)  # single batch random sampling for debugging purposes

    def __getitem__(self, item):
        return self.data_loaders_dict[item]

    def __setitem__(self, key, value):
        self.data_loaders_dict[key] = value

    def __contains__(self, item):
        return item in self.data_loaders_dict

    def __iter__(self):
        return zip_first_cycle(*[self[name] for name in self.data_loaders_loop])

    def __call__(self, shuffle=False, indices=None):
        if indices is not None and shuffle:
            raise ValueError('indices is mutually exclusive with shuffle')
        if indices is None:
            if shuffle:
                sampler = RandomSampler(self.gene_dataset)
            else:
                sampler = SequentialSampler(self.gene_dataset)
        else:
            sampler = SubsetRandomSampler(indices)
        return DataLoader(self.gene_dataset, sampler=sampler, **self.data_loaders_kwargs)

    @staticmethod
    def raw_data(*data_loaders):
        """
        From a list of data_loaders, return the numpy array suitable for classification with sklearn
        :param data_loaders: a sequence of data_loaders arguments
        :return: a sequence of (data, labels) for each data_loader
        """

        def get_data_labels(data_loader):
            if hasattr(data_loader, 'sampler') and hasattr(data_loader.sampler, 'indices'):
                data = data_loader.dataset.X[data_loader.sampler.indices]
                labels = data_loader.dataset.labels[data_loader.sampler.indices]
            else:
                data = data_loader.dataset.X
                labels = data_loader.dataset.labels
            return data, labels.ravel()

        return [get_data_labels(data_loader) for data_loader in data_loaders]


class TrainTestDataLoaders(DataLoaders):
    to_monitor = ['train', 'test']
    data_loaders_loop = ['train']

    def __init__(self, gene_dataset, train_size=0.1, test_size=None, seed=0, **data_loaders_kwargs):
        """
        :param train_size: float, int, or None (default is 0.1)
        :param test_size: float, int, or None (default is None)
        """
        super(TrainTestDataLoaders, self).__init__(gene_dataset, **data_loaders_kwargs)

        n = len(self.gene_dataset)
        n_train, n_test = _validate_shuffle_split(n, test_size, train_size)
        np.random.seed(seed=seed)
        permutation = np.random.permutation(n)
        indices_test = permutation[:n_test]
        indices_train = permutation[n_test:(n_test + n_train)]

        data_loader_train = self(indices=indices_train)
        data_loader_test = self(indices=indices_test)

        self.data_loaders_dict.update({
            'train': data_loader_train,
            'test': data_loader_test
        })


class SemiSupervisedDataLoaders(DataLoaders):
    to_monitor = ['labelled', 'unlabelled']

    def __init__(self, gene_dataset, n_labelled_samples_per_class=50, seed=0, **data_loaders_kwargs):
        """
        :param n_labelled_samples_per_class: number of labelled samples per class
        """
        super(SemiSupervisedDataLoaders, self).__init__(gene_dataset, **data_loaders_kwargs)
        indices_labelled, indices_unlabelled = (
            get_indices(gene_dataset.labels, [n_labelled_samples_per_class] * gene_dataset.n_labels, seed=seed)
        )
        data_loader_all = self(shuffle=True)
        data_loader_labelled = self(indices=indices_labelled)
        data_loader_unlabelled = self(indices=indices_unlabelled)

        self.data_loaders_dict.update({
            'all': data_loader_all,
            'labelled': data_loader_labelled,
            'unlabelled': data_loader_unlabelled,
        })


class JointSemiSupervisedDataLoaders(SemiSupervisedDataLoaders):
    data_loaders_loop = ['all', 'labelled']


class AlternateSemiSupervisedDataLoaders(SemiSupervisedDataLoaders):
    data_loaders_loop = ['all']

    def classifier_data_loaders(self):
        data_loaders = DataLoaders(gene_dataset=self.gene_dataset, **self.data_loaders_kwargs)
        data_loaders.data_loaders_dict.update({
            'train': self['labelled'],
            'test': self['unlabelled']
        })
        data_loaders.to_monitor = ['train']
        return data_loaders
