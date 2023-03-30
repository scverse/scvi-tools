from docrep import DocstringProcessor

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

param_shuffle = """\
shuffle
    Whether or not to shuffle the data before splitting. Ignored if `train_indices` is
    not `None`."""

param_pin_memory = """\
pin_memory
    Whether to copy tensors into device-pinned memory before returning them. Passed
    into :class:`~scvi.data.AnnDataLoader`. Will take priority over the value of
    `scvi.settings.dl_pin_memory_gpu_training`."""


data_splitting_dsp = DocstringProcessor(
    param_n_obs=param_n_obs,
    param_train_size=param_train_size,
    param_validation_size=param_validation_size,
    param_train_indices=param_train_indices,
    param_validation_indices=param_validation_indices,
    param_shuffle=param_shuffle,
)
