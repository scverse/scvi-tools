"""Shared docstrings."""

doc_differential_expression = """\
adata
    AnnData object with equivalent structure to initial AnnData.
    If None, defaults to the AnnData object used to initialize the model.
groupby
    The key of the observations grouping to consider.
group1
    Subset of groups, e.g. [`'g1'`, `'g2'`, `'g3'`], to which comparison
    shall be restricted, or all groups in `groupby` (default).
group2
    If `None`, compare each group in `group1` to the union of the rest of the groups
    in `groupby`. If a group identifier, compare with respect to this group.
idx1
    Boolean mask or indices for `group1`. `idx1` and `idx2` can be used as an alternative
    to the AnnData keys. If `idx1` is not `None`, this option overrides `group1`
    and `group2`.
idx2
    Boolean mask or indices for `group2`. By default, includes all cells not specified in
    `idx1`.
mode
    Method for differential expression. See user guide for full explanation.
delta
    specific case of region inducing differential expression.
    In this case, we suppose that :math:`R \setminus [-\delta, \delta]` does not induce differential expression
    (change model default case).
batch_size
    Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
all_stats
    Concatenate count statistics (e.g., mean expression group 1) to DE results.
batch_correction
    Whether to correct for batch effects in DE inference.
batchid1
    Subset of categories from `batch_key` registered in :func:`~scvi.data.setup_anndata`,
    e.g. [`'batch1'`, `'batch2'`, `'batch3'`], for `group1`. Only used if `batch_correction` is `True`, and
    by default all categories are used.
batchid2
    Same as `batchid1` for group2. `batchid2` must either have null intersection with `batchid1`,
    or be exactly equal to `batchid1`. When the two sets are exactly equal, cells are compared by
    decoding on the same batch. When sets have null intersection, cells from `group1` and `group2`
    are decoded on each group in `group1` and `group2`, respectively.
fdr_target
    Tag features as DE based on posterior expected false discovery rate.
"""
