from docrep import DocstringProcessor

de_adata = """\
adata
    AnnData object with equivalent structure to initial AnnData.
    If None, defaults to the AnnData object used to initialize the model."""
de_groupby = """\
groupby
    The key of the observations grouping to consider."""
de_group1 = """\
group1
    Subset of groups, e.g. [`'g1'`, `'g2'`, `'g3'`], to which comparison
    shall be restricted, or all groups in `groupby` (default)."""
de_group2 = """\
group2
    If `None`, compare each group in `group1` to the union of the rest of the groups
    in `groupby`. If a group identifier, compare with respect to this group."""
de_idx1 = """\
idx1
    `idx1` and `idx2` can be used as an alternative to the AnnData keys.
    Custom identifier for `group1` that can be of three sorts: (1) a boolean mask,
    (2) indices, or (3) a string. If it is a string, then it will query indices that
    verifies conditions on `adata.obs`, as described in :meth:`pandas.DataFrame.query`
    If `idx1` is not `None`, this option overrides `group1`
    and `group2`."""
de_idx2 = """\
idx2
    Custom identifier for `group2` that has the same
    properties as `idx1`.
    By default, includes all cells not specified in
    `idx1`."""
de_mode = """\
mode
    Method for differential expression. See user guide for full explanation."""
de_delta = """\
delta
    specific case of region inducing differential expression.
    In this case, we suppose that :math:`R \\setminus [-\\delta, \\delta]` does not induce differential expression
    (change model default case)."""
de_batch_size = """\
batch_size
    Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`."""
de_all_stats = """\
all_stats
    Concatenate count statistics (e.g., mean expression group 1) to DE results."""
de_batch_correction = """\
batch_correction
    Whether to correct for batch effects in DE inference."""
de_batchid1 = """\
batchid1
    Subset of categories from `batch_key` registered in ``setup_anndata``,
    e.g. [`'batch1'`, `'batch2'`, `'batch3'`], for `group1`. Only used if `batch_correction` is `True`, and
    by default all categories are used."""
de_batchid2 = """\
batchid2
    Same as `batchid1` for group2. `batchid2` must either have null intersection with `batchid1`,
    or be exactly equal to `batchid1`. When the two sets are exactly equal, cells are compared by
    decoding on the same batch. When sets have null intersection, cells from `group1` and `group2`
    are decoded on each group in `group1` and `group2`, respectively."""
de_fdr_target = """\
fdr_target
    Tag features as DE based on posterior expected false discovery rate."""
de_silent = """\
silent
    If True, disables the progress bar. Default: False."""
de_importance_sampling = """\
importance_sampling
    Whether to use importance sampling to compute normalized gene expression."""
de_fn_kwargs = """\
fn_kwargs
    Additional kwargs for the normalized gene expression estimation.
    Only applies if `importance_sampling` is True."""

de_dsp = DocstringProcessor(
    de_adata=de_adata,
    de_groupby=de_groupby,
    de_group1=de_group1,
    de_group2=de_group2,
    de_idx1=de_idx1,
    de_idx2=de_idx2,
    de_mode=de_mode,
    de_delta=de_delta,
    de_batch_size=de_batch_size,
    de_all_stats=de_all_stats,
    de_batch_correction=de_batch_correction,
    de_batchid1=de_batchid1,
    de_batchid2=de_batchid2,
    de_fdr_target=de_fdr_target,
    de_silent=de_silent,
    de_importance_sampling=de_importance_sampling,
    de_fn_kwargs=de_fn_kwargs,
)


summary = """\
Sets up the :class:`~anndata.AnnData` object for this model.

A mapping will be created between data fields used by this model to their respective locations in adata.
None of the data in adata are modified. Only adds fields to adata"""

summary_mdata = """\
Sets up the :class:`~mudata.MuData` object for this model.

A mapping will be created between data fields used by this model to their respective locations in adata.
None of the data in adata are modified. Only adds fields to adata"""

param_mdata = """\
mdata
    MuData object. Rows represent cells, columns represent features."""

param_adata = """\
adata
    AnnData object. Rows represent cells, columns represent features."""

param_batch_key = """\
batch_key
    key in `adata.obs` for batch information. Categories will automatically be converted into integer
    categories and saved to `adata.obs['_scvi_batch']`. If `None`, assigns the same batch to all the data."""

param_labels_key = """\
labels_key
    key in `adata.obs` for label information. Categories will automatically be converted into integer
    categories and saved to `adata.obs['_scvi_labels']`. If `None`, assigns the same label to all the data."""

param_size_factor_key = """\
size_factor_key
    key in `adata.obs` for size factor information. Instead of using library size as a size factor, the provided
    size factor column will be used as offset in the mean of the likelihood. Assumed to be on linear scale."""

param_layer = """\
layer
    if not `None`, uses this as the key in `adata.layers` for raw count data."""

param_cat_cov_keys = """\
categorical_covariate_keys
    keys in `adata.obs` that correspond to categorical data.
    These covariates can be added in addition to the batch covariate and are also treated as nuisance factors
    (i.e., the model tries to minimize their effects on the latent space). Thus, these should not be used for
    biologically-relevant factors that you do _not_ want to correct for."""

param_cont_cov_keys = """\
continuous_covariate_keys
    keys in `adata.obs` that correspond to continuous data.
    These covariates can be added in addition to the batch covariate and are also treated as nuisance factors
    (i.e., the model tries to minimize their effects on the latent space). Thus, these should not be used for
    biologically-relevant factors that you do _not_ want to correct for."""

param_unlabeled_category = """\
unlabeled_category
    value in `adata.obs[labels_key]` that indicates unlabeled observations."""

param_modalities = """\
modalities
    Dictionary mapping parameters to modalities."""

param_copy = """\
copy
    if `True`, a copy of adata is returned."""

returns = """\
None. Adds the following fields:

.uns['_scvi']
    `scvi` setup dictionary
.obs['_scvi_labels']
    labels encoded as integers
.obs['_scvi_batch']
    batch encoded as integers"""


setup_anndata_dsp = DocstringProcessor(
    summary=summary,
    summary_mdata=summary_mdata,
    param_mdata=param_mdata,
    param_adata=param_adata,
    param_batch_key=param_batch_key,
    param_labels_key=param_labels_key,
    param_layer=param_layer,
    param_cat_cov_keys=param_cat_cov_keys,
    param_cont_cov_keys=param_cont_cov_keys,
    param_size_factor_key=param_size_factor_key,
    param_unlabeled_category=param_unlabeled_category,
    param_modalities=param_modalities,
    param_copy=param_copy,
    returns=returns,
)


param_accelerator = """\
accelerator
    Supports passing different accelerator types `("cpu", "gpu", "tpu", "ipu", "hpu",
    "mps, "auto")` as well as custom accelerator instances."""

param_devices = """\
devices
    The devices to use. Can be set to a non-negative index (`int` or `str`), a sequence
    of device indices (`list` or comma-separated `str`), the value `-1` to indicate all
    available devices, or `"auto"` for automatic selection based on the chosen
    `accelerator`. If set to `"auto"` and `accelerator` is not determined to be `"cpu"`,
    then `devices` will be set to the first available device."""

param_device = """\
device
    The device to use. Can be set to a non-negative index (`int` or `str`) or `"auto"`
    for automatic selection based on the chosen accelerator. If set to `"auto"` and
    `accelerator` is not determined to be `"cpu"`, then `device` will be set to the
    first available device."""

param_return_device = """\
return_device
    Returns the first or only device as determined by `accelerator` and `devices`.
    Depending on the value, will either return a PyTorch device (`"torch"`), a Jax
    device (`"jax"`), or neither (`None`)."""

param_validate_single_device = """\
validate_single_device
    Validates that `devices` is set to a single device if `device!="auto"` and throws
    an error if not."""


devices_dsp = DocstringProcessor(
    param_accelerator=param_accelerator,
    param_devices=param_devices,
    param_device=param_device,
    param_return_device=param_return_device,
    param_validate_single_device=param_validate_single_device,
)
