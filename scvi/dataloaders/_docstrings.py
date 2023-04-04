from docrep import DocstringProcessor

from scvi.utils._docstrings import param_accelerator, param_device

summary = """\
Creates train, validation, and test split dataloaders."""

param_adata_manager = """\
adata_manager
    :class:`~scvi.data.AnnDataManager` object that has been created via a model's
    `setup_anndata` or `setup_mudata` method."""

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

param_pin_memory = """\
pin_memory
    Whether to copy tensors into page-locked memory before returning them, which enables
    faster data transfer to CUDA-enabled GPUs. Passed into
    :class:`~scvi.data.AnnDataLoader`. Equivalent to setting
    `scvi.settings.dl_pin_memory_gpu_training`."""

param_n_samples_per_label = """\
n_samples_per_label
    Number of subsamples for each label class to sample per epoch."""

param_dataloader_kwargs = """\
train_dataloader_kwargs
    Keyword arguments passed into the dataloaders, an instance of
    :class:`~scvi.data.AnnDataLoader` if there is no labeled data or
    :class:`~scvi.data.SemiSupervisedDataLoader` if there is labeled data."""

param_device_backed = """\
device_backed
    If `True`, the data will be loaded into the device specified by `accelerator` and
    `device` prior to training. This can speed up training, but requires more memory."""


dataloaders_dsp = DocstringProcessor(
    summary=summary,
    param_adata_manager=param_adata_manager,
    param_n_obs=param_n_obs,
    param_train_size=param_train_size,
    param_validation_size=param_validation_size,
    param_all_indices=param_all_indices,
    param_train_indices=param_train_indices,
    param_validation_indices=param_validation_indices,
    param_shuffle=param_datasplitter_shuffle,
    param_pin_memory=param_pin_memory,
    param_n_samples_per_label=param_n_samples_per_label,
    param_dataloader_kwargs=param_dataloader_kwargs,
    param_accelerator=param_accelerator,
    param_device=param_device,
)

dataloader_summary = """\
Loads :class:`~torch.Tensor`s from an :class:`~anndata.AnnData` dataset."""

param_indices = """\
indices
    The indices of the observations in the dataset to include. If `None`, all
    observations will be included."""

param_indices_list = """\
indices_list
    A list where each element is a list of the indices of the observations in the dataset
    to include."""

param_dataloader_shuffle = """\
shuffle
    Whether or not to shuffle the data for sampling."""

param_batch_size = """\
batch_size
    The size of the minibatch to load in each iteration."""

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

param_drop_last = """\
drop_last
    If `True`, the last incomplete minibatch will be dropped if the number of
    observations is not divisible by the batch size."""

param_iter_ndarray = """\
iter_ndarray
    If `True`, the dataloader will iterate over the data as a numpy array."""

param_kwargs = """\
**kwargs
    Keyword arguments passed into :class:`~torch.utils.data.DataLoader`."""


dataloader_dsp = DocstringProcessor(
    summary=dataloader_summary,
    param_adata_manager=param_adata_manager,
    param_shuffle=param_dataloader_shuffle,
    param_indices=param_indices,
    param_indices_list=param_indices_list,
    param_n_samples_per_label=param_n_samples_per_label,
    param_batch_size=param_batch_size,
    param_data_and_attributes=param_data_and_attributes,
    param_drop_last=param_drop_last,
    param_iter_ndarray=param_iter_ndarray,
    param_accelerator=param_accelerator,
    param_device=param_device,
    param_device_backed=param_device_backed,
    param_kwargs=param_kwargs,
)
