import copy
from itertools import cycle

import numpy as np
import torch
from sklearn.model_selection._split import _validate_shuffle_split
from torch.utils.data import DataLoader
from torch.utils.data.sampler import SubsetRandomSampler, SequentialSampler, RandomSampler


class SequentialSubsetSampler(SubsetRandomSampler):
    def __init__(self, indices):
        self.indices = np.sort(indices)

    def __iter__(self):
        return iter(self.indices)


class DataLoaderWrapper(DataLoader):
    def __init__(self, dataset, use_cuda=True, **kwargs):
        self.kwargs = kwargs
        self.use_cuda = use_cuda and torch.cuda.is_available()
        super(DataLoaderWrapper, self).__init__(dataset, **kwargs)

    def to_cuda(self, tensors):
        return [t.cuda(async=self.use_cuda) if self.use_cuda else t for t in tensors]

    def sequential(self, batch_size=128):
        kwargs = copy.copy(self.kwargs)
        kwargs['batch_size'] = batch_size
        if hasattr(self, 'sampler') and hasattr(self.sampler, 'indices'):
            kwargs['sampler'] = SequentialSubsetSampler(indices=self.sampler.indices)
        else:
            kwargs['sampler'] = SequentialSampler(self.dataset)
        return DataLoaderWrapper(self.dataset, use_cuda=self.use_cuda, **kwargs)

    def __iter__(self):
        return map(self.to_cuda, super(DataLoaderWrapper, self).__iter__())


class DataLoaders:
    to_monitor = []
    loop = []

    def __init__(self, gene_dataset, use_cuda=True, **kwargs):
        """
        :param gene_dataset: a GeneExpressionDataset instance
        :param kwargs: any additional keyword arguments to pass to the data_loaders at .init
        """
        self.gene_dataset = gene_dataset
        self.use_cuda = use_cuda
        self.kwargs = {
            "batch_size": 128,
            "pin_memory": use_cuda,
            "collate_fn": self.gene_dataset.collate_fn
        }
        self.kwargs.update(kwargs)
        self.dict = {'sequential': self()}
        self.infinite_random_sampler = cycle(self(shuffle=True))  # convenient for debugging - see .sample()

    def sample(self):
        return next(self.infinite_random_sampler)  # single batch random sampling for debugging purposes

    def __getitem__(self, item):
        return self.dict[item]

    def __setitem__(self, key, value):
        self.dict[key] = value

    def __contains__(self, item):
        return item in self.dict

    def __iter__(self):
        data_loaders_loop = [self[name] for name in self.loop]
        return zip(data_loaders_loop[0], *[cycle(data_loader) for data_loader in data_loaders_loop[1:]])

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
        return DataLoaderWrapper(self.gene_dataset, use_cuda=self.use_cuda, sampler=sampler,
                                 **self.kwargs)

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
    loop = ['train']

    def __init__(self, gene_dataset, train_size=0.1, test_size=None, seed=0, **kwargs):
        """
        :param train_size: float, int, or None (default is 0.1)
        :param test_size: float, int, or None (default is None)
        """
        super(TrainTestDataLoaders, self).__init__(gene_dataset, **kwargs)

        n = len(self.gene_dataset)
        n_train, n_test = _validate_shuffle_split(n, test_size, train_size)
        np.random.seed(seed=seed)
        permutation = np.random.permutation(n)
        indices_test = permutation[:n_test]
        indices_train = permutation[n_test:(n_test + n_train)]

        data_loader_train = self(indices=indices_train)
        data_loader_test = self(indices=indices_test)

        self.dict.update({
            'train': data_loader_train,
            'test': data_loader_test
        })


class SemiSupervisedDataLoaders(DataLoaders):
    to_monitor = ['labelled', 'unlabelled']

    def __init__(self, gene_dataset, n_labelled_samples_per_class=50, seed=0, use_cuda=True, **kwargs):
        """
        :param n_labelled_samples_per_class: number of labelled samples per class
        """
        super(SemiSupervisedDataLoaders, self).__init__(gene_dataset, use_cuda=use_cuda, **kwargs)

        n_labelled_samples_per_class_array = [n_labelled_samples_per_class] * gene_dataset.n_labels
        labels = np.array(gene_dataset.labels).ravel()
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
        indices_labelled = permutation_idx[indices[:total_labelled]]
        indices_unlabelled = permutation_idx[indices[total_labelled:]]

        data_loader_all = self(shuffle=True)
        data_loader_labelled = self(indices=indices_labelled)
        data_loader_unlabelled = self(indices=indices_unlabelled)

        self.dict.update({
            'all': data_loader_all,
            'labelled': data_loader_labelled,
            'unlabelled': data_loader_unlabelled,
        })


class JointSemiSupervisedDataLoaders(SemiSupervisedDataLoaders):
    loop = ['all', 'labelled']


class AlternateSemiSupervisedDataLoaders(SemiSupervisedDataLoaders):
    loop = ['all']

    def classifier_data_loaders(self):
        data_loaders = DataLoaders(gene_dataset=self.gene_dataset, use_cuda=self.use_cuda, **self.kwargs)
        data_loaders.loop = ['train']
        data_loaders.dict.update({
            'train': self['labelled'],
            'test': self['unlabelled']
        })
        data_loaders.to_monitor = ['train']
        return data_loaders
