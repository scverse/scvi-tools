from docrep import DocstringProcessor

from scvi.utils._docstrings import param_accelerator, param_device

dataset_summary = """\
A :class:`~torch.utils.data.Dataset` that loads tensors from :class:`~anndata.AnnData`."""

param_adata_manager = """\
adata_manager
    :class:`~scvi.data.AnnDataManager` object that has been created via a model's
    `setup_anndata` or `setup_mudata` method."""

param_data_and_attributes = """\
data_and_attributes
    Used if not all fields registered with `setup_adata` or `setup_mudata` should be
    loaded in training minibatches. One of the following:

    * A dictionary with keys corresponding to a subset of the keys in the data registry
    (`adata_manager.data_registry.keys`) to be loaded in each minibatch, and values
    corresponding to thedesired loading data types.

    * A list of keys corresponding to a subset of the keys in the data registry
    (`adata_manager.data_registry.keys`) to be loaded in each minibatch.

    * `None`: All registered data will be loaded in each minibatch."""

param_device_backed = """\
device_backed
    If `True`, the data will be loaded into the device specified by `accelerator` and
    `device` prior to training. This can speed up training, but it is only recommended
    if the entire data fits in device memory."""

dataset_dsp = DocstringProcessor(
    summary=dataset_summary,
    param_adata_manager=param_adata_manager,
    param_getitem_tensors=param_data_and_attributes,
    param_accelerator=param_accelerator,
    param_device=param_device,
    param_device_backed=param_device_backed,
)

dataloader_summary = """\
A :class:`~torch.utils.data.DataLoader` that loads from :class:`~anndata.AnnData`."""

param_dataset = """\
dataset
    The :class:`~torch.utils.data.Dataset` to load from. Cannot be set if `adata_manager`
    is set. Must be set if `adata_manager` is not set."""

param_indices = """\
indices
    The indices of the observations in the dataset to include. If `None`, all
    observations will be included."""

param_indices_list = """\
indices_list
    A list where each element is a list of the indices of the observations in the dataset
    to include."""

param_sampler = """\
sampler
    Defines the strategy to draw samples from the dataset. Can be any `Iterable` with
    `__len__` implemented. If specified, will override the default batch sampler."""

param_dataloader_shuffle = """\
shuffle
    If `True`, :class:`~torch.utils.data.RandomSampler` will be used to sample
    observations for minibatches. If `False`, :class:`~torch.utils.data.SequentialSampler`
    will be used instead."""

param_batch_size = """\
batch_size
    The size of the minibatch to load in each iteration."""

param_n_samples_per_label = """\
n_samples_per_label
    The size of the subset to sample from each set of unique labels per training epoch.
    If `None`, the dataloader will sample from all labeled observations in each epoch."""

param_drop_last = """\
drop_last
    If `True`, the last incomplete minibatch will be dropped if the number of
    observations is not divisible by the batch size."""

param_iter_ndarray = """\
iter_ndarray
    If `True`, the dataloader will iterate over the data as a numpy array."""

param_seed = """\
seed
    Random seed for the resampling generator."""

param_pin_memory = """\
pin_memory
    Whether to copy tensors into page-locked memory before returning them, which enables
    faster data transfer to CUDA-enabled GPUs. Passed into
    :class:`~scvi.data.AnnDataLoader`. Equivalent to setting
    `scvi.settings.dl_pin_memory_gpu_training`."""

param_dataloader_kwargs = """\
**kwargs
    Keyword arguments passed into :class:`~torch.utils.data.DataLoader`."""

dataloader_dsp = DocstringProcessor(
    summary=dataloader_summary,
    param_adata_manager=param_adata_manager,
    param_dataset=param_dataset,
    param_indices=param_indices,
    param_indices_list=param_indices_list,
    param_sampler=param_sampler,
    param_shuffle=param_dataloader_shuffle,
    param_n_samples_per_label=param_n_samples_per_label,
    param_batch_size=param_batch_size,
    param_data_and_attributes=param_data_and_attributes,
    param_drop_last=param_drop_last,
    param_iter_ndarray=param_iter_ndarray,
    param_accelerator=param_accelerator,
    param_device=param_device,
    param_device_backed=param_device_backed,
    param_pin_memory=param_pin_memory,
    param_kwargs=param_dataloader_kwargs,
)

datasplitter_summary = """\
Creates train, validation, and test split dataloaders."""

param_n_obs = """\
n_obs
    The number of observations in the dataset. Must be `None` if `all_indices` is not
    `None`."""

param_train_size = """\
train_size
    The fraction of the dataset to assign to the training set. Must be `None` if
    `train_indices` is not `None`."""

param_validation_size = """\
validation_size
    The fraction of the dataset to assign to the validation set. If `None`, the
    validation set will correspond to the remaining data. Must be `None` if `train_size`
    is `None`. If `train_size + validation_size < 1`, the remaining observations
    will be assigned to the test set."""

param_all_indices = """\
all_indices
    The indices of all the observations in the dataset. Must be `None` if `n_obs` is not
    `None`."""

param_train_indices = """\
train_indices
    The indices of the observations belonging to the training set. Must be `None` if
    `train_size` is not `None`."""

param_validation_indices = """\
validation_indices
    The indices of the observations belonging to the validation set. Must be `None` if
    `train_indices is `None`. If `None`, the validation set will correspond to the
    remaining data. If the union of the two sets of indices is not the full set of
    indices, the remaining indices will be assigned to the test set."""

param_datasplitter_shuffle = """\
shuffle
    Whether or not to shuffle the data before splitting. Ignored if `train_indices` is
    not `None`."""

param_datasplitter_kwargs = """\
train_dataloader_kwargs
    Keyword arguments passed into the dataloaders, an instance of
    :class:`~scvi.data.AnnDataLoader` if there is no labeled data or
    :class:`~scvi.data.SemiSupervisedDataLoader` if there is labeled data."""

datasplitter_dsp = DocstringProcessor(
    summary=datasplitter_summary,
    param_adata_manager=param_adata_manager,
    param_n_obs=param_n_obs,
    param_train_size=param_train_size,
    param_validation_size=param_validation_size,
    param_all_indices=param_all_indices,
    param_train_indices=param_train_indices,
    param_validation_indices=param_validation_indices,
    param_shuffle=param_datasplitter_shuffle,
    param_n_samples_per_label=param_n_samples_per_label,
    param_pin_memory=param_pin_memory,
    param_kwargs=param_datasplitter_kwargs,
)
