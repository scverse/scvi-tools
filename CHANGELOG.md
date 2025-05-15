# Release notes

Starting from version 0.20.1, this format is based on [Keep a Changelog], and this project adheres
to [Semantic Versioning]. Full commit history is available in the
[commit logs](https://github.com/scverse/scvi-tools/commits/).

## Version 1.3

### 1.3.1 (2025-05-15)

#### Added

- Add {class}`scvi.external.METHYLANVI` for modeling methylation labelled data from single-cell
    bisulfite sequencing (scBS-seq) {pr}`3066`.
- Add supervised module class {class}`scvi.module.base.SupervisedModuleClass`. {pr}`3237`.
- Add get normalized function model property for any generative model {pr}`3238` and changed
    get_accessibility_estimates to get_normalized_accessibility, where needed.
- Add {class}`scvi.external.TOTALANVI`. {pr}`3259`.
- Add Custom Dataloaders registry support, {pr}`2932`.
- Add support for using Census and LaminAI custom dataloaders for {class}`scvi.model.SCVI`
    and {class}`scvi.model.SCANVI`, {pr}`2932`.
- Add Early stopping KL warmup steps. {pr}`3262`.
- Add Minification option to {class}`~scvi.model.LinearSCVI` {pr}`3294`.
- Update Read the docs tutorials index page with interactive filterable options {pr}`3276`.

#### Fixed

- Add consideration for missing monitor set during early stopping. {pr}`3226`.
- Fix bug in SysVI get_normalized_expression function. {pr}`3255`.
- Add support for IntegratedGradients for multimodal models. {pr}`3264`.
- Fix bug in resolVI get_normalized expression function. {pr}`3308`.
- Fix bug in resolVI gene-assay dispersion. {pr}`3308`.

#### Changed

- Updated Scvi-Tools AWS hub to Weizmann instead of Berkeley. {pr}`3246`.
- Updated resolVI to use rapids-singlecell. {pr}`3308`.

#### Removed

- Removed Jax version constraint for mrVI training. {pr}`3309`.

### 1.3.0 (2025-02-28)

#### Added

- Add {class}`scvi.external.Decipher` for dimensionality reduction and interpretable
    representation learning in single-cell RNA sequencing data {pr}`3015`, {pr}`3091`.
- Add multiGPU support for {class}`~scvi.model.SCVI`, {class}`~scvi.model.SCANVI`,
    {class}`~scvi.model.CondSCVI` and {class}`~scvi.model.LinearSCVI`, {class}`~scvi.model.TOTALVI`,
    {class}`~scvi.model.MULTIVI` and {class}`~scvi.model.PEAKVI`. {pr}`3125`.
- Add an exception callback to {class}`scvi.train._callbacks.SaveCheckpoint` in order to save
    optimal model during training, in case of failure because of Nan's in gradients. {pr}`3159`.
- Add {meth}`~scvi.model.SCVI.get_normalized_expression` for models: {class}`~scvi.model.PEAKVI`,
    {class}`~scvi.external.POISSONVI`, {class}`~scvi.model.CondSCVI`, {class}`~scvi.model.AUTOZI`,
    {class}`~scvi.external.CellAssign` and {class}`~scvi.external.GIMVI`. {pr}`3121`.
- Add {class}`scvi.external.RESOLVI` for bias correction in single-cell resolved spatial
    transcriptomics {pr}`3144`.
- Add semisupervised training mixin class
    {class}`scvi.model.base.SemisupervisedTrainingMixin`. {pr}`3164`.
- Add scib-metrics support for {class}`scvi.autotune.AutotuneExperiment` and
    {class}`scvi.train._callbacks.ScibCallback` for autotune for scib metrics {pr}`3168`.
- Add Support of dask arrays in AnnTorchDataset. {pr}`3193`.
- Add a common use cases section in the docs user guide. {pr}`3200`.
- Add {class}`scvi.external.SysVI` for cycle consistency loss and VampPrior {pr}`3195`.

#### Fixed

- Fixed bug in distributed {class}`scvi.dataloaders.ConcatDataLoader` {pr}`3053`.
- Fixed bug when loading Pyro-based models and scArches support for Pyro {pr}`3138`
- Fixed disable vmap in {class}`scvi.external.MRVI` for large sample sizes to avoid
    out-of-memory errors. Store distance matrices as numpy array in xarray to reduce
    memory usage {pr}`3146`.
- Fixed {class}`scvi.external.MRVI` MixtureSameFamily log probability calculation {pr}`3189`.

#### Changed

- Updated the CI workflow with multiGPU tests {pr}`3053`.
- Set `mode="change"` as default DE method. Compute positive and negative LFC separately
    by default (`test_mode="three"`). Corrected computation of pseudocounts and make if
    default to add a pseudocounts for genes not expressed (`pseudocount=None`). According to
    Eq. 10 of Boyeau _et al_, _PNAS_ 2023 {pr}`2826`

#### Removed

## Version 1.2

### 1.2.2 (2024-12-31)

#### Added

- Add MuData Minification option to {class}`~scvi.model.TOTALVI` {pr}`3061`.
- Add Support for MPS usage in mac {pr}`3100`.
- Add support for torch.compile before train (EXPERIMENTAL) {pr}`2931`.
- Add support for Numpy 2.0 {pr}`2842`.
- Changed scvi-hub ModelCard and add criticism metrics to the card {pr}`3078`.
- MuData support for {class}`~scvi.model.MULTIVI` via the method
    {meth}`~scvi.model.MULTIVI.setup_mudata` {pr}`3038`.

#### Fixed

- Fixed batch_size pop to get in {class}`scvi.dataloaders.DataSplitter` {pr}`3128`.

#### Changed

- Updated the CI workflow with internet, private and optional tests {pr}`3082`.
- Changed loompy stored files to anndata {pr}`2842`.
- Address AnnData >= 0.11 deprecation warning for {class}`anndata.experimental` by replacing
    instances to {class}`anndata.abc` and {class}`anndata.io` {pr}`3085`.

#### Removed

- Removed the support for loompy and local mde function {pr}`2842`.

### 1.2.1 (2024-12-04)

#### Added

- Added adaptive handling for last training minibatch of 1-2 cells in case of
    `datasplitter_kwargs={"drop_last": False}` and `train_size = None` by moving them into
    validation set, if available. {pr}`3036`.
- Add `batch_key` and `labels_key` to {meth}`scvi.external.SCAR.setup_anndata`. {pr}`3045`.
- Implemented variance of ZINB distribution. {pr}`3044`.
- Support for minified mode while retaining counts to skip the encoder.
- New Trainingplan argument `update_only_decoder` to use stored latent codes and skip training of
    the encoder.
- Refactored code for minified models. {pr}`2883`.
- Add {class}`scvi.external.METHYLVI` for modeling methylation data from single-cell
    bisulfite sequencing (scBS-seq) experiments {pr}`2834`.

#### Fixed

- Breaking Change: Fix `get_outlier_cell_sample_pairs` function in {class}`scvi.external.MRVI`
    to correctly compute the maxmimum log-density across in-sample cells rather than the
    aggregated posterior log-density {pr}`3007`.
- Fix references to `scvi.external` in {meth}`scvi.external.SCAR.setup_anndata`.
- Fix gimVI to append mini batches first into CPU during get_imputed and get_latent operations {pr}`3058`.

#### Changed

#### Removed

### 1.2.0 (2024-09-26)

#### Added

- Add support for Python 3.12 {pr}`2966`.
- Add support for categorial covariates in scArches in {class}`scvi.model.base.ArchesMixin` {pr}`2936`.
- Add assertion error in cellAssign for checking duplicates in celltype markers {pr}`2951`.
- Add {meth}`scvi.external.POISSONVI.get_region_factors` {pr}`2940`.
- {attr}`scvi.settings.dl_persistent_workers` allows using persistent workers in
    {class}`scvi.dataloaders.AnnDataLoader` {pr}`2924`.
- Add option for using external indexes in data splitting classes that are under `scvi.dataloaders`
    by passing `external_indexing=list[train_idx,valid_idx,test_idx]` as well as in all models
    available {pr}`2902`.
- Add warning if creating data splits in `scvi.dataloaders` that create last batch with less than 3
    cells {pr}`2916`.
- Add new experimental functional API for hyperparameter tuning with
    {func}`scvi.autotune.run_autotune` and {class}`scvi.autotune.AutotuneExperiment` to replace
    {class}`scvi.autotune.ModelTuner`, {class}`scvi.autotune.TunerManager`, and
    {class}`scvi.autotune.TuneAnalysis` {pr}`2561`.
- Add experimental class {class}`scvi.nn.Embedding` implementing methods for extending embeddings
    {pr}`2574`.
- Add experimental support for representing batches with continuously-valued embeddings by passing
    in `batch_representation="embedding"` to {class}`scvi.model.SCVI` {pr}`2576`.
- Add experimental mixin classes {class}`scvi.model.base.EmbeddingMixin` and
    {class}`scvi.module.base.EmbeddingModuleMixin` {pr}`2576`.
- Add option to generate synthetic spatial coordinates in {func}`scvi.data.synthetic_iid` with
    argument `generate_coordinates` {pr}`2603`.
- Add experimental support for using custom {class}`lightning.pytorch.core.LightningDataModule`s
    in {func}`scvi.autotune.run_autotune` {pr}`2605`.
- Add {class}`scvi.external.VELOVI` for RNA velocity estimation using variational inference
    {pr}`2611`.
- Add `unsigned` argument to {meth}`scvi.hub.HubModel.pull_from_s3` to allow for unsigned
    downloads of models from AWS S3 {pr}`2615`.
- Add support for `batch_key` in {meth}`scvi.model.CondSCVI.setup_anndata` {pr}`2626`.
- Add support for {meth}`scvi.model.base.RNASeqMixin` in {class}`scvi.model.CondSCVI` {pr}`2915`.
- Add `load_best_on_end` argument to {class}`scvi.train.SaveCheckpoint` to load the best model
    state at the end of training {pr}`2672`.
- Add experimental class {class}`scvi.distributions.BetaBinomial` implementing the Beta-Binomial
    distribution with mean-dispersion parameterization for modeling scBS-seq methylation data
    {pr}`2692`.
- Add support for custom dataloaders in {class}`scvi.model.base.VAEMixin` methods by specifying
    the `dataloader` argument {pr}`2748`.
- Add option to use a normal distribution in the generative model of {class}`scvi.model.SCVI` by
    passing in `gene_likelihood="normal"` {pr}`2780`.
- Add {class}`scvi.external.MRVI` for modeling sample-level heterogeneity in single-cell RNA-seq
    data {pr}`2756`.
- Add support for reference mapping with {class}`mudata.MuData` models to
    {class}`scvi.model.base.ArchesMixin` {pr}`2578`.
- Add argument `return_mean` to {meth}`scvi.model.base.VAEMixin.get_reconstruction_error`
    and {meth}`scvi.model.base.VAEMixin.get_elbo` to allow computation
    without averaging across cells {pr}`2362`.
- Add support for setting `weights="importance"` in
    {meth}`scvi.model.SCANVI.differential_expression` {pr}`2362`.

#### Changed

- Deprecate {func}`scvi.data.cellxgene`, to be removed in v1.3. Please directly use the
    [cellxgene-census](https://chanzuckerberg.github.io/cellxgene-census/) instead {pr}`2542`.
- Deprecate {func}`scvi.nn.one_hot`, to be removed in v1.3. Please directly use the
    `one_hot` function in PyTorch instead {pr}`2608`.
- Deprecate {class}`scvi.train.SaveBestState`, to be removed in v1.3. Please use
    {class}`scvi.train.SaveCheckpoint` instead {pr}`2673`.
- Deprecate `save_best` argument in {meth}`scvi.model.PEAKVI.train` and
    {meth}`scvi.model.MULTIVI.train`, to be removed in v1.3. Please pass in `enable_checkpointing`
    or specify a custom checkpointing procedure with {class}`scvi.train.SaveCheckpoint` instead
    {pr}`2673`.
- Move {func}`scvi.model.base._utils._load_legacy_saved_files` to
    {func}`scvi.model.base._save_load._load_legacy_saved_files` {pr}`2731`.
- Move {func}`scvi.model.base._utils._load_saved_files` to
    {func}`scvi.model.base._save_load._load_saved_files` {pr}`2731`.
- Move {func}`scvi.model.base._utils._initialize_model` to
    {func}`scvi.model.base._save_load._initialize_model` {pr}`2731`.
- Move {func}`scvi.model.base._utils._validate_var_names` to
    {func}`scvi.model.base._save_load._validate_var_names` {pr}`2731`.
- Move {func}`scvi.model.base._utils._prepare_obs` to
    {func}`scvi.model.base._de_core._prepare_obs` {pr}`2731`.
- Move {func}`scvi.model.base._utils._de_core` to
    {func}`scvi.model.base._de_core._de_core` {pr}`2731`.
- Move {func}`scvi.model.base._utils._fdr_de_prediction` to
    {func}`scvi.model.base._de_core_._fdr_de_prediction` {pr}`2731`.
- {func}`scvi.data.synthetic_iid` now generates unique variable names for protein and
    accessibility data {pr}`2739`.
- The `data_module` argument in {meth}`scvi.model.base.UnsupervisedTrainingMixin.train` has been
    renamed to `datamodule` for consistency {pr}`2749`.
- Change the default saving method of variable names for {class}`mudata.MuData` based models
    (_e.g._ {class}`scvi.model.TOTALVI`) to a dictionary of per-mod variable names instead of a
    concatenated array of all variable names. Users may replicate the previous behavior by
    passing in `legacy_mudata_format=True` to {meth}`scvi.model.base.BaseModelClass.save`
    {pr}`2769`.
- Changed internal activation function in {class}`scvi.nn.DecoderTOTALVI` to Softplus to
    increase numerical stability. This is the new default for new models. Previously trained models
    will be loaded with exponential activation function {pr}`2913`.

#### Fixed

- Fix logging of accuracy for cases with 1 sample per class in scANVI {pr}`2938`.
- Disable adversarial classifier if training with a single batch.
    Previously this raised a None error {pr}`2914`.
- {meth}`~scvi.model.SCVI.get_normalized_expression` fixed for Poisson distribution and
    Negative Binomial with latent_library_size {pr}`2915`.
- Fix {meth}`scvi.module.VAE.marginal_ll` when `n_mc_samples_per_pass=1` {pr}`2362`.
- {meth}`scvi.module.VAE.marginal_ll` when `n_mc_samples_per_pass=1` {pr}`2362`.
- Enable option to drop_last minibatch during training by `datasplitter_kwargs={"drop_last": True}`
    {pr}`2926`.
- Fix JAX to be deterministic on CUDA when seed is manually set {pr}`2923`.

#### Removed

- Remove {class}`scvi.autotune.ModelTuner`, {class}`scvi.autotune.TunerManager`, and
    {class}`scvi.autotune.TuneAnalysis` in favor of new experimental functional API with
    {func}`scvi.autotune.run_autotune` and {class}`scvi.autotune.AutotuneExperiment` {pr}`2561`.
- Remove `feed_labels` argument and corresponding code paths in {meth}`scvi.module.SCANVAE.loss`
    {pr}`2644`.
- Remove {class}`scvi.train._callbacks.MetricsCallback` and argument `additional_val_metrics` in
    {class}`scvi.train.Trainer` {pr}`2646`.

## Version 1.1

### 1.1.6 (2024-08-19)

#### Fixed

- Breaking change: In `scvi.autotune._manager` we changed the parameter in RunConfig from
    `local_dir` to `storage_path` see issue `2908` {pr}`2689`.

### 1.1.5 (2024-06-30)

### 1.1.4 (2024-06-30)

#### Added

- Add argument `return_logits` to {meth}`scvi.external.SOLO.predict` that allows returning logits
    instead of probabilities when passing in `soft=True` to replicate the buggy behavior previous
    to v1.1.3 {pr}`2870`.

### 1.1.3 (2024-06-26)

#### Fixed

- Breaking change: Fix {meth}`scvi.external.SOLO.predict` to correctly return probabiities
    instead of logits when passing in `soft=True` (the default option) {pr}`2689`.
- Breaking change: Fix {class}`scvi.dataloaders.SemiSupervisedDataSplitter` to properly sample
    unlabeled observations without replacement {pr}`2816`.

### 1.1.2 (2024-03-01)

#### Changed

- Address AnnData >= 0.10 deprecation warning for {func}`anndata.read` by replacing instances with
    {func}`anndata.read_h5ad` {pr}`2531`.
- Address AnnData >= 0.10 deprecation warning for {class}`anndata._core.sparse_dataset.SparseDataset`
    by replacing instances with {class}`anndata.experimental.CSCDataset` and
    {class}`anndata.experimental.CSRDataset` {pr}`2531`.

### 1.1.1 (2024-02-19)

#### Fixed

- Correctly apply non-default user parameters in {class}`scvi.external.POISSONVI` {pr}`2522`.

### 1.1.0 (2024-02-13)

#### Added

- Add {class}`scvi.external.ContrastiveVI` for contrastiveVI {pr}`2242`.
- Add {class}`scvi.dataloaders.BatchDistributedSampler` for distributed training {pr}`2102`.
- Add `additional_val_metrics` argument to {class}`scvi.train.Trainer`, allowing to specify
    additional metrics to compute and log during the validation loop using
    {class}`scvi.train._callbacks.MetricsCallback` {pr}`2136`.
- Expose `accelerator` and `device` arguments in {meth}`scvi.hub.HubModel.load_model` `pr`{2166}.
- Add `load_sparse_tensor` argument in {class}`scvi.data.AnnTorchDataset` for directly loading
    SciPy CSR and CSC data structures to their PyTorch counterparts, leading to faster data loading
    depending on the sparsity of the data {pr}`2158`.
- Add per-group LFC information to
    {meth}`scvi.criticism.PosteriorPredictiveCheck.differential_expression`. `metrics["diff_exp"]`
    is now a dictionary where `summary` stores the summary dataframe, and `lfc_per_model_per_group`
    stores the per-group LFC {pr}`2173`.
- Expose {meth}`torch.save` keyword arguments in {class}`scvi.model.base.BaseModelClass.save`
    and {class}`scvi.external.GIMVI.save` {pr}`2200`.
- Add `model_kwargs` and `train_kwargs` arguments to {meth}`scvi.autotune.ModelTuner.fit`
    {pr}`2203`.
- Add `datasplitter_kwargs` to model `train` methods {pr}`2204`.
- Add `use_posterior_mean` argument to {meth}`scvi.model.SCANVI.predict` for stochastic prediction
    of celltype labels {pr}`2224`.
- Add support for Python 3.10+ type annotations in {class}`scvi.autotune.ModelTuner` {pr}`2239`.
- Add option to log device statistics in {meth}`scvi.autotune.ModelTuner.fit` with argument
    `monitor_device_stats` {pr}`2260`.
- Add option to pass in a random seed to {meth}`scvi.autotune.ModelTuner.fit` with argument `seed`
    {pr}`2260`.
- Automatically log the learning rate when `reduce_lr_on_plateau=True` in training plans
    {pr}`2280`.
- Add {class}`scvi.external.POISSONVI` to model scATAC-seq fragment counts with a Poisson
    distribution {pr}`2249`
- {class}`scvi.train.SemiSupervisedTrainingPlan` now logs the classifier calibration error
    {pr}`2299`.
- Passing `enable_checkpointing=True` into `train` methods is now compatible with our model saves.
    Additional options can be specified by initializing with {class}`scvi.train.SaveCheckpoint`
    {pr}`2317`.
- {attr}`scvi.settings.dl_num_workers` is now correctly applied as the default `num_workers` in
    {class}`scvi.dataloaders.AnnDataLoader` {pr}`2322`.
- Passing in `indices` to {class}`scvi.criticism.PosteriorPredictiveCheck` allows for running
    metrics on a subset of the data {pr}`2361`.
- Add `seed` argument to {func}`scvi.model.utils.mde` for reproducibility {pr}`2373`.
- Add {meth}`scvi.hub.HubModel.save` and {meth}`scvi.hub.HubMetadata.save` {pr}`2382`.
- Add support for Optax 0.1.8 by renaming instances of {func}`optax.additive_weight_decay` to
    {func}`optax.add_weight_decay` {pr}`2396`.
- Add support for hosting {class}`scvi.hub.HubModel` on AWS S3 via
    {meth}`scvi.hub.HubModel.pull_from_s3` and {meth}`scvi.hub.HubModel.push_to_s3` {pr}`2378`.
- Add clearer error message for {func}`scvi.data.poisson_gene_selection` when input data does not
    contain raw counts {pr}`2422`.
- Add API for using custom dataloaders with {class}`scvi.model.SCVI` by making `adata` argument
    optional on initialization and adding optional argument `data_module` to
    {meth}`scvi.model.base.UnsupervisedTrainingMixin.train` {pr}`2467`.
- Add support for Ray 2.8-2.9 in {class}`scvi.autotune.ModelTuner` {pr}`2478`.

#### Fixed

- Fix bug where `n_hidden` was not being passed into {class}`scvi.nn.Encoder` in
    {class}`scvi.model.AmortizedLDA` {pr}`2229`
- Fix bug in {class}`scvi.module.SCANVAE` where classifier probabilities were interpreted as
    logits. This is backwards compatible as loading older models will use the old code path
    {pr}`2301`.
- Fix bug in {class}`scvi.external.GIMVI` where `batch_size` was not properly used in inference
    methods {pr}`2366`.
- Fix error message formatting in {meth}`scvi.data.fields.LayerField.transfer_field` {pr}`2368`.
- Fix ambiguous error raised in {meth}`scvi.distributions.NegativeBinomial.log_prob` and
    {meth}`scvi.distributions.ZeroInflatedNegativeBinomial.log_prob` when `scale` not passed in
    and value not in support {pr}`2395`.
- Fix initialization of {class}`scvi.distributions.NegativeBinomial` and
    {class}`scvi.distributions.ZeroInflatedNegativeBinomial` when `validate_args=True` and
    optional parameters not passed in {pr}`2395`.
- Fix error when re-initializing {class}`scvi.external.GIMVI` with the same datasets {pr}`2446`.

#### Changed

- Replace `sparse` with `sparse_format` argument in {meth}`scvi.data.synthetic_iid` for increased
    flexibility over dataset format {pr}`2163`.
- Revalidate `devices` when automatically switching from MPS to CPU accelerator in
    {func}`scvi.model._utils.parse_device_args` {pr}`2247`.
- Refactor {class}`scvi.data.AnnTorchDataset`, now loads continuous data as {class}`numpy.float32`
    and categorical data as {class}`numpy.int64` by default {pr}`2250`.
- Support fractional GPU usage in {class}`scvi.autotune.ModelTuner` `pr`{2252}.
- Tensorboard is now the default logger in {class}`scvi.autotune.ModelTuner` `pr`{2260}.
- Match `momentum` and `epsilon` in {class}`scvi.module.JaxVAE` to the default values in PyTorch
    {pr}`2309`.
- Change {class}`scvi.train.SemiSupervisedTrainingPlan` and
    {class}`scvi.train.ClassifierTrainingPlan` accuracy and F1 score
    computations to use `"micro"` reduction rather than `"macro"` {pr}`2339`.
- Internal refactoring of {meth}`scvi.module.VAE.sample` and
    {meth}`scvi.model.base.RNASeqMixin.posterior_predictive_sample` {pr}`2377`.
- Change `xarray` and `sparse` from mandatory to optional dependencies {pr}`2480`.
- Use {class}`anndata.experimental.CSCDataset` and {class}`anndata.experimental.CSRDataset`
    instead of the deprecated {class}`anndata._core.sparse_dataset.SparseDataset` for type checks
    {pr}`2485`.
- Make `use_observed_lib_size` argument adjustable in {class}`scvi.module.LDVAE` `pr`{2494}.

#### Removed

- Remove deprecated `use_gpu` argument in favor of PyTorch Lightning arguments `accelerator` and
    `devices` {pr}`2114`.
- Remove deprecated `scvi._compat.Literal` class {pr}`2115`.
- Remove chex dependency {pr}`2482`.

## Version 1.0

### 1.0.4 (2023-10-13)

### Added

- Add support for AnnData 0.10.0 {pr}`2271`.

### 1.0.3 (2023-08-13)

#### Changed

- Disable the default selection of MPS when `accelerator="auto"` in Lightning {pr}`2167`.
- Change JAX models to use `dict` instead of {class}`flax.core.FrozenDict` according
    to the Flax migration guide <https://github.com/google/flax/discussions/3191> {pr}`2222`.

#### Fixed

- Fix bug in {class}`scvi.model.base.PyroSviTrainMixin` where `training_plan`
    argument is ignored {pr}`2162`.
- Fix missing docstring for `unlabeled_category` in
    {class}`scvi.model.SCANVI.setup_anndata` and reorder arguments {pr}`2189`.
- Fix Pandas 2.0 unpickling error in {meth}`scvi.model.base.BaseModelClas.convert_legacy_save`
    by switching to {func}`pandas.read_pickle` for the setup dictionary {pr}`2212`.

### 1.0.2 (2023-07-05)

#### Fixed

- Fix link to Scanpy preprocessing in introduction tutorial {pr}`2154`.
- Fix link to Ray Tune search API in autotune tutorial {pr}`2154`.

### 1.0.1 (2023-07-04)

#### Added

- Add support for Python 3.11 {pr}`1977`.

#### Changed

- Upper bound Chex dependency to 0.1.8 due to NumPy installation conflicts {pr}`2132`.

### 1.0.0 (2023-06-02)

#### Added

- Add {class}`scvi.criticism.PosteriorPredictiveCheck` for model evaluation {pr}`2058`.
- Add {func}`scvi.data.reads_to_fragments` for scATAC data {pr}`1946`
- Add default `stacklevel` for `warnings` in `scvi.settings` {pr}`1971`.
- Add scBasset motif injection procedure {pr}`2010`.
- Add importance sampling based differential expression procedure {pr}`1872`.
- Raise clearer error when initializing {class}`scvi.external.SOLO` from {class}`scvi.model.SCVI`
    with extra categorical or continuous covariates {pr}`2027`.
- Add option to generate {class}`mudata.MuData` in {meth}`scvi.data.synthetic_iid` {pr}`2028`.
- Add option for disabling shuffling prior to splitting data in
    {class}`scvi.dataloaders.DataSplitter` {pr}`2037`.
- Add {meth}`scvi.data.AnnDataManager.create_torch_dataset` and expose custom sampler ability
    {pr}`2036`.
- Log training loss through Lightning's progress bar {pr}`2043`.
- Filter Jax undetected GPU warnings {pr}`2044`.
- Raise warning if MPS backend is selected for PyTorch,
    see <https://github.com/pytorch/pytorch/issues/77764> {pr}`2045`.
- Add `deregister_manager` function to {class}`scvi.model.base.BaseModelClass`, allowing to clear
    {class}`scvi.data.AnnDataManager` instances from memory {pr}`2060`.
- Add option to use a linear classifier in {class}`scvi.model.SCANVI` {pr}`2063`.
- Add lower bound 0.12.1 for Numpyro dependency {pr}`2078`.
- Add new section in scBasset tutorial for motif scoring {pr}`2079`.

#### Fixed

- Fix creation of minified adata by copying original uns dict {pr}`2000`. This issue arises with
    anndata>=0.9.0.
- Fix {class}`scvi.model.TOTALVI` {class}`scvi.model.MULTIVI` handling of missing protein values
    {pr}`2009`.
- Fix bug in {meth}`scvi.distributions.NegativeBinomialMixture.sample` where `theta` and `mu`
    arguments were switched around {pr}`2024`.
- Fix bug in {meth}`scvi.dataloaders.SemiSupervisedDataLoader.resample_labels` where the labeled
    dataloader was not being reinitialized on subsample {pr}`2032`.
- Fix typo in {class}`scvi.model.JaxSCVI` example snippet {pr}`2075`.

#### Changed

- Use sphinx book theme for documentation {pr}`1673`.
- {meth}`scvi.model.base.RNASeqMixin.posterior_predictive_sample` now outputs 3-d
    {class}`sparse.GCXS` matrices {pr}`1902`.
- Add an option to specify `dropout_ratio` in {meth}`scvi.data.synthetic_iid` {pr}`1920`.
- Update to lightning 2.0 {pr}`1961`
- Hyperopt is new default searcher for tuner {pr}`1961`
- {class}`scvi.train.AdversarialTrainingPlan` no longer encodes data twice during a training step,
    instead uses same latent for both optimizers {pr}`1961`, {pr}`1980`
- Switch back to using sphinx autodoc typehints {pr}`1970`.
- Disable default seed, run `scvi.settings.seed` after import for reproducibility {pr}`1976`.
- Deprecate `use_gpu` in favor of PyTorch Lightning arguments `accelerator` and `devices`, to be
    removed in v1.1 {pr}`1978`.
- Docs organization {pr}`1983`.
- Validate training data and code URLs for {class}`scvi.hub.HubMetadata` and
    {class}`scvi.hub.HubModelCardHelper` {pr}`1985`.
- Keyword arguments for encoders and decoders can now be passed in from the model level {pr}`1986`.
- Expose `local_dir` as a public property in {class}`scvi.hub.HubModel` {pr}`1994`.
- Use {func}`anndata.concat` internally inside {meth}`scvi.external.SOLO.from_scvi_model` {pr}`2013`.
- {class}`scvi.train.SemiSupervisedTrainingPlan` and {class}`scvi.train.ClassifierTrainingPlan`
    now log accuracy, F1 score, and AUROC metrics {pr}`2023`.
- Switch to cellxgene census for backend for cellxgene data function {pr}`2030`.
- Change default `max_cells` and `truncation` in
    {meth}`scvi.model.base.RNASeqMixin._get_importance_weights` {pr}`2064`.
- Refactor heuristic for default `max_epochs` as a separate function
    {meth}`scvi.model._utils.get_max_epochs_heuristic` {pr}`2083`.

#### Removed

- Remove ability to set up ST data in {class}`~scvi.external.SpatialStereoscope.from_rna_model`,
    which was deprecated. ST data should be set up using
    {class}`~scvi.external.SpatialStereoscope.setup_anndata` {pr}`1949`.
- Remove custom reusable doc decorator which was used for de docs {pr}`1970`.
- Remove `drop_last` as an integer from {class}`~scvi.dataloaders.AnnDataLoader`, add typing and
    code cleanup {pr}`1975`.
- Remove seqfish and seqfish plus datasets {pr}`2017`.
- Remove support for Python 3.8 (NEP 29) {pr}`2021`.

## Version 0.20

### 0.20.3 (2023-03-21)

#### Fixed

- Fix totalVI differential expression when integer sequential protein names are automatically used
    {pr}`1951`.
- Fix peakVI scArches test case {pr}`1962`.

#### Changed

- Allow passing in `map_location` into {meth}`~scvi.hub.HubMetadata.from_dir` and
    {meth}`~scvi.hub.HubModelCardHelper.from_dir` and set default to `"cpu"` {pr}`1960`.
- Updated tutorials {pr}`1966`.

### 0.20.2 (2023-03-10)

#### Fixed

- Fix `return_dist` docstring of {meth}`scvi.model.base.VAEMixin.get_latent_representation`
    {pr}`1932`.
- Fix hyperlink to pymde docs {pr}`1944`

#### Changed

- Use ruff for fixing and linting {pr}`1921`, {pr}`1941`.
- Use sphinx autodoc instead of sphinx-autodoc-typehints {pr}`1941`.
- Remove .flake8 and .prospector files {pr}`1923`.
- Log individual loss terms in {meth}`scvi.module.MULTIVAE.loss` {pr}`1936`.
- Setting up ST data in {class}`~scvi.external.SpatialStereoscope.from_rna_model` is deprecated.
    ST data should be set up using {class}`~scvi.external.SpatialStereoscope.setup_anndata`
    {pr}`1803`.

### 0.20.1 (2023-02-21)

#### Fixed

- Fixed computation of ELBO during training plan logging when using global kl terms. {pr}`1895`
- Fixed usage of {class}`scvi.train.SaveBestState` callback, which affected
    {class}`scvi.model.PEAKVI` training. If using {class}`~scvi.model.PEAKVI`, please upgrade.
    {pr}`1913`
- Fixed original seed for jax-based models to work with jax 0.4.4. {pr}`1907`, {pr}`1909`

### New in 0.20.0 (2023-02-01)

#### Major changes

- Model hyperparameter tuning is available through {class}`~scvi.autotune.ModelTuner` (beta)
    {pr}`1785`,{pr}`1802`,{pr}`1831`.
- Pre-trained models can now be uploaded to and downloaded from [Hugging Face models] using the
    {mod}`~scvi.hub` module {pr}`1779`,{pr}`1812`,{pr}`1828`,{pr}`1841`, {pr}`1851`,{pr}`1862`.
- {class}`~anndata.AnnData` `.var` and `.varm` attributes can now be registered through new fields
    in {mod}`~scvi.data.fields` {pr}`1830`,{pr}`1839`.
- {class}`~scvi.external.SCBASSET`, a reimplementation of the [original scBasset model], is
    available for representation learning of scATAC-seq data (experimental) {pr}`1839`,{pr}`1844`,
    {pr}`1867`,{pr}`1874`,{pr}`1882`.
- {class}`~scvi.train.LowLevelPyroTrainingPlan` and {class}`~scvi.model.base.PyroModelGuideWarmup`
    added to allow the use of vanilla PyTorch optimization on Pyro models {pr}`1845`,{pr}`1847`.
- Add {meth}`scvi.data.cellxgene` function to download cellxgene datasets {pr}`1880`.

#### Minor changes

- Latent mode support changed so that user data is no longer edited in-place {pr}`1756`.
- Minimum supported PyTorch Lightning version is now 1.9 {pr}`1795`,{pr}`1833`,{pr}`1863`.
- Minimum supported Python version is now 3.8 {pr}`1819`.
- [Poetry] removed in favor of [Hatch] for builds and publishing {pr}`1823`.
- `setup_anndata` docstrings fixed, `setup_mudata` docstrings added {pr}`1834`,{pr}`1837`.
- {meth}`~scvi.data.add_dna_sequence` adds DNA sequences to {class}`~anndata.AnnData` objects using
    [genomepy] {pr}`1839`,{pr}`1842`.
- Update tutorial formatting with pre-commit {pr}`1850`
- Expose `accelerators` and `devices` arguments in {class}`~scvi.train.Trainer` {pr}`1864`.
- Development in GitHub Codespaces is now supported {pr}`1836`.

#### Breaking changes

- {class}`~scvi.module.base.LossRecorder` has been removed in favor of
    {class}`~scvi.module.base.LossOutput` {pr}`1869`.

#### Bug Fixes

- {class}`~scvi.train.JaxTrainingPlan` now correctly updates `global_step` through PyTorch
    Lightning by using a dummy optimizer. {pr}`1791`.
- CUDA compatibility issue fixed in {meth}`~scvi.distributions.ZeroInflatedNegativeBinomial.sample`
    {pr}`1813`.
- Device-backed {class}`~scvi.dataloaders.AnnTorchDataset` fixed to work with sparse data {pr}`1824`.
- Fix bug {meth}`~scvi.model.base._log_likelihood.compute_reconstruction_error` causing the first
    batch to be ignored, see more details in {issue}`1854` {pr}`1857`.

#### Contributors

- {ghuser}`adamgayoso`
- {ghuser}`eroell`
- {ghuser}`gokceneraslan`
- {ghuser}`macwiatrak`
- {ghuser}`martinkim0`
- {ghuser}`saroudant`
- {ghuser}`vitkl`
- {ghuser}`watiss`

## Version 0.19

### New in 0.19.0 (2022-10-31)

#### Major Changes

- {class}`~scvi.train.TrainingPlan` allows custom PyTorch optimizers [#1747].
- Improvements to {class}`~scvi.train.JaxTrainingPlan` [#1747] [#1749].
- {class}`~scvi.module.base.LossRecorder` is deprecated. Please substitute with
    {class}`~scvi.module.base.LossOutput` [#1749]
- All training plans require keyword args after the first positional argument [#1749]
- {class}`~scvi.module.base.JaxBaseModuleClass` absorbed features from the `JaxModuleWrapper`,
    rendering the `JaxModuleWrapper` obsolote, so it was removed. [#1751]
- Add {class}`scvi.external.Tangram` and {class}`scvi.external.tangram.TangramMapper` that
    implement Tangram for mapping scRNA-seq data to spatial data [#1743].

#### Minor changes

- Remove confusing warning about kl warmup, log kl weight instead [#1773]

#### Breaking changes

- {class}`~scvi.module.base.LossRecorder` no longer allows access to dictionaries of values if
    provided during initialization [#1749].
- `JaxModuleWrapper` removed. [#1751]

#### Bug Fixes

- Fix `n_proteins` usage in {class}`~scvi.model.MULTIVI` [#1737].
- Remove unused param in {class}`~scvi.model.MULTIVI` [#1741].
- Fix random seed handling for Jax models [#1751].

#### Contributors

- [@watiss]
- [@adamgayoso]
- [@martinkim0]
- [@marianogabitto]

## Version 0.18

### New in 0.18.0 (2022-10-12)

#### Major Changes

- Add latent mode support in {class}`~scvi.model.SCVI` [#1672]. This allows for loading a model
    using latent representations only (i.e. without the full counts). Not only does this speed up
    inference by using the cached latent distribution parameters (thus skipping the encoding step),
    but this also helps in scenarios where the full counts are not available but cached latent
    parameters are. We provide utility functions and methods to dynamically convert a model to
    latent mode.
- Added {class}`~scvi.external.SCAR` as an external model for ambient RNA removal [#1683].

#### Minor changes

- Faster inference in PyTorch with `torch.inference_mode` [#1695].
- Upgrade to Lightning 1.6 [#1719].
- Update CI workflow to separate static code checking from pytest [#1710].
- Add Python 3.10 to CI workflow [#1711].
- Add {meth}`~scvi.data.AnnDataManager.register_new_fields` [#1689].
- Use sphinxcontrib-bibtex for references [#1731].
- {meth}`~scvi.model.base.VAEMixin.get_latent_representation`: more explicit and better docstring
    [#1732].
- Replace custom attrdict with {class}`~ml_collections` implementation [#1696].

#### Breaking changes

- Add weight support to {class}`~scvi.model.MULTIVI` [#1697]. Old models can't be loaded anymore.

#### Bug Fixes

- Fix links for breast cancer and mouse datasets [#1709].
- fix quick start notebooks not showing [#1733].

#### Contributors

- [@watiss]
- [@adamgayoso]
- [@martinkim0]
- [@ricomnl]
- [@marianogabitto]

## Version 0.17

### New in 0.17.4 (2021-09-20)

#### Changes

- Support for PyTorch Lightning 1.7 [#1622].
- Allow `flax` to use any mutable states used by a model generically with
    {class}`~scvi.module.base.TrainStateWithState` [#1665], [#1700].
- Update publication links in `README` [#1667].
- Docs now include floating window cross references with `hoverxref`, external links with
    `linkcode`, and `grid` [#1678].

#### Bug Fixes

- Fix `get_likelihood_parameters()` failure when `gene_likelihood != "zinb"` in
    {class}`~scvi.model.base.RNASeqMixin` [#1618].
- Fix exception logic when not using the observed library size in {class}`~scvi.module.VAE`
    initialization [#1660].
- Replace instances of `super().__init__()` with an argument in `super()`, causing `autoreload`
    extension to throw errors [#1671].
- Change cell2location tutorial causing docs build to fail [#1674].
- Replace instances of `max_epochs` as `int`s for new PyTorch Lightning [#1686].
- Catch case when `torch.backends.mps` is not implemented [#1692].
- Fix Poisson sampling in {meth}`~scvi.module.VAE.sample` [#1702].

#### Contributors

- [@adamgayoso]
- [@watiss]
- [@mkarikom]
- [@tommycelsius]
- [@ricomnl]

### New in 0.17.3 (2022-08-26)

#### Changes

- Pin sphinx_gallery to fix tutorial cards on docs [#1657]
- Use latest tutorials in release [#1657]

#### Contributors

- [@watiss]
- [@adamgayoso]

### New in 0.17.2 (2022-08-26)

#### Changes

- Move `training` argument in {class}`~scvi.module.JaxVAE` constructor to a keyword argument into
    the call method. This simplifies the {class}`~scvi.module.base.JaxModuleWrapper` logic and
    avoids the reinstantiation of {class}`~scvi.module.JaxVAE` during evaluation [#1580].
- Add a static method on the BaseModelClass to return the AnnDataManger's full registry [#1617].
- Clarify docstrings for continuous and categorical covariate keys [#1637].
- Remove poetry lock, use newer build system [#1645].

#### Bug Fixes

- Fix CellAssign to accept extra categorical covariates [#1629].
- Fix an issue where `max_epochs` is never determined heuristically for totalvi, instead it would
    always default to 400 [#1639].

#### Breaking Changes

- Fix an issue where `max_epochs` is never determined heuristically for totalvi, instead it would
    always default to 400 [#1639].

#### Contributors

- [@watiss]
- [@RK900]
- [@adamgayoso]
- [@jjhong922]

### New in 0.17.1 (2022-07-14)

Make sure notebooks are up to date for real this time :).

#### Contributors

- [@jjhong922]
- [@adamgayoso]

### New in 0.17.0 (2022-07-14)

#### Major Changes

- Experimental MuData support for {class}`~scvi.model.TOTALVI` via the method
    {meth}`~scvi.model.TOTALVI.setup_mudata`. For several of the existing `AnnDataField` classes,
    there is now a MuData counterpart with an additional `mod_key` argument used to indicate the
    modality where the data lives (e.g. {class}`~scvi.data.fields.LayerField` to
    {class}`~scvi.data.fields.MuDataLayerField`). These modified classes are simply wrapped
    versions of the original `AnnDataField` code via the new
    {class}`scvi.data.fields.MuDataWrapper` method [#1474].

- Modification of the {meth}`~scvi.module.VAE.generative` method's outputs to return prior and
    likelihood properties as {class}`~torch.distributions.distribution.Distribution` objects.
    Concerned modules are {class}`~scvi.module.AmortizedLDAPyroModule`, {class}`AutoZIVAE`,
    {class}`~scvi.module.MULTIVAE`, {class}`~scvi.module.PEAKVAE`, {class}`~scvi.module.TOTALVAE`,
    {class}`~scvi.module.SCANVAE`, {class}`~scvi.module.VAE`, and {class}`~scvi.module.VAEC`. This
    allows facilitating the manipulation of these distributions for model training and inference
    [#1356].

- Major changes to Jax support for scvi-tools models to generalize beyond
    {class}`~scvi.model.JaxSCVI`. Support for Jax remains experimental and is subject to breaking
    changes:

    - Consistent module interface for Flax modules (Jax-backed) via
        {class}`~scvi.module.base.JaxModuleWrapper`, such that they are compatible with the
        existing {class}`~scvi.model.base.BaseModelClass` [#1506].
    - {class}`~scvi.train.JaxTrainingPlan` now leverages Pytorch Lightning to factor out
        Jax-specific training loop implementation [#1506].
    - Enable basic device management in Jax-backed modules [#1585].

#### Minor changes

- Add {meth}`~scvi.module.base.PyroBaseModuleClass.on_load` callback which is called on
    {meth}`~scvi.model.base.BaseModuleClass.load` prior to loading the module state dict [#1542].
- Refactor metrics code and use {class}`~torchmetrics.MetricCollection` to update metrics in bulk
    [#1529].
- Add `max_kl_weight` and `min_kl_weight` to {class}`~scvi.train.TrainingPlan` [#1595].
- Add a warning to {class}`~scvi.model.base.UnsupervisedTrainingMixin` that is raised if
    `max_kl_weight` is not reached during training [#1595].

#### Breaking changes

- Any methods relying on the output of `inference` and `generative` from existing scvi-tools models
    (e.g. {class}`~scvi.model.SCVI`, {class}`~scvi.model.SCANVI`) will need to be modified to
    accept `torch.Distribution` objects rather than tensors for each parameter (e.g. `px_m`,
    `px_v`) [#1356].
- The signature of {meth}`~scvi.train.TrainingPlan.compute_and_log_metrics` has changed to support
    the use of {class}`~torchmetrics.MetricCollection`. The typical modification required will look
    like changing `self.compute_and_log_metrics(scvi_loss, self.elbo_train)` to
    `self.compute_and_log_metrics(scvi_loss, self.train_metrics, "train")`. The same is necessary
    for validation metrics except with `self.val_metrics` and the mode `"validation"` [#1529].

#### Bug Fixes

- Fix issue with {meth}`~scvi.model.SCVI.get_normalized_expression` with multiple samples and
    additional continuous covariates. This bug originated from {meth}`~scvi.module.VAE.generative`
    failing to match the dimensions of the continuous covariates with the input when `n_samples>1`
    in {meth}`~scvi.module.VAE.inference` in multiple module classes [#1548].
- Add support for padding layers in {meth}`~scvi.model.SCVI.prepare_query_anndata` which is
    necessary to run {meth}`~scvi.model.SCVI.load_query_data` for a model setup with a layer
    instead of X [#1575].

#### Contributors

- [@jjhong922]
- [@adamgayoso]
- [@PierreBoyeau]
- [@RK900]
- [@FlorianBarkmann]

## Version 0.16

### New in 0.16.4 (2022-06-14)

Note: When applying any model using the {class}`~scvi.train.AdversarialTrainingPlan` (e.g.
{class}`~scvi.model.TOTALVI`, {class}`~scvi.model.MULTIVI`), you should make sure to use v0.16.4
instead of v0.16.3 or v0.16.2. This release fixes a critical bug in the training plan.

#### Changes

#### Breaking changes

#### Bug Fixes

- Fix critical issue in {class}`~scvi.train.AdversarialTrainingPlan` where `kl_weight` was
    overwritten to 0 at each step ([#1566]). Users should avoid using v0.16.2 and v0.16.3 which
    both include this bug.

#### Contributors

- [@jjhong922]
- [@adamgayoso]

### New in 0.16.3 (2022-06-04)

#### Changes

- Removes sphinx max version and removes jinja dependency ([#1555]).

#### Breaking changes

#### Bug Fixes

- Upper bounds protobuf due to pytorch lightning incompatibilities ([#1556]). Note that [#1556]
    has unique changes as PyTorch Lightning >=1.6.4 adds the upper bound in their requirements.

#### Contributors

- [@jjhong922]
- [@adamgayoso]

### New in 0.16.2 (2022-05-10)

#### Changes

#### Breaking changes

#### Bug Fixes

- Raise appropriate error when `backup_url` is not provided and file is missing on
    {meth}`~scvi.model.base.BaseModelClass.load` ([#1527]).
- Pipe `loss_kwargs` properly in {class}`~scvi.train.AdversarialTrainingPlan`, and fix incorrectly
    piped kwargs in {class}`~scvi.model.TOTALVI` and {class}`~scvi.model.MULTIVI` ([#1532]).

#### Contributors

- [@jjhong922]
- [@adamgayoso]

### New in 0.16.1 (2022-04-22)

#### Changes

- Update scArches Pancreas tutorial, DestVI tutorial ([#1520]).

#### Breaking changes

- {class}`~scvi.dataloaders.SemiSupervisedDataLoader` and
    {class}`~scvi.dataloaders.SemiSupervisedDataSplitter` no longer take `unlabeled_category` as an
    initial argument. Instead, the `unlabeled_category` is fetched from the labels state registry,
    assuming that the {class}`~scvi.data.AnnDataManager` object is registered with a
    {class}`~scvi.data.fields.LabelsWithUnlabeledObsField` ([#1515]).

#### Bug Fixes

- Bug fixed in {class}`~scvi.model.SCANVI` where `self._labeled_indices` was being improperly set
    ([#1515]).
- Fix issue where {class}`~scvi.model.SCANVI.load_query_data` would not properly add an obs column
    with the unlabeled category when the `labels_key` was not present in the query data.
- Disable extension of categories for labels in {class}`~scvi.model.SCANVI.load_query_data`
    ([#1519]).
- Fix an issue with {meth}`~scvi.model.SCANVI.prepare_query_data` to ensure it does nothing when
    genes are completely matched ([#1520]).

#### Contributors

- [@jjhong922]
- [@adamgayoso]

### New in 0.16.0 (2022-04-12)

This release features a refactor of {class}`~scvi.model.DestVI` ([#1457]):

1. Bug fix in cell type amortization, which leads to on par performance of cell type amortization
    `V_encoder` with free parameter for cell type proportions `V`.
1. Bug fix in library size in {class}`~scvi.model.CondSCVI`, that lead to downstream dependency
    between sum over cell type proportions `v_ind` and library size `library` in
    {class}`~scvi.model.DestVI`.
1. `neg_log_likelihood_prior` is not computed anymore on random subset of single cells but cell
    type specific subclustering using cluster variance `var_vprior`, cluster mean `mean_vprior` and
    cluster mixture proportion `mp_vprior` for computation. This leads to more stable results and
    faster computation time. Setting `vamp_prior_p` in {func}`~scvi.model.DestVI.from_rna_model` to
    the expected resolution is critical in this algorithm.
1. The new default is to also use dropout `dropout` during the decoder of
    {class}`~scvi.model.CondSCVI` and subsequently `dropout_decoder` in
    {class}`~scvi.model.DestVI`, we found this to be beneficial after bug fixes listed above.
1. We changed the weighting of the loss on the variances of beta and the prior of eta.

:::{note}
Due to bug fixes listed above this version of {class}`~scvi.model.DestVI` is not backwards
compatible. Despite instability in training in the outdated version, we were able to reproduce
results generated with this code. We therefore do not strictly encourage to rerun old experiments.
:::

We published a new tutorial. This new tutorial incorporates a new utility package
[destvi_utils](https://github.com/YosefLab/destvi_utils) that generates exploratory plots of the
results of {class}`~scvi.model.DestVI`. We refer to the manual of this package for further
documentation.

#### Changes

- Docs changes (installation [#1498], {class}`~scvi.model.DestVI` user guide [#1501] and [#1508],
    dark mode code cells [#1499]).
- Add `backup_url` to the {meth}`~scvi.model.base.BaseModelClass.load` method of each model class,
    enabling automatic downloading of model save file ([#1505]).

#### Breaking changes

- Support for loading legacy loading is removed from {meth}`~scvi.model.base.BaseModelClass.load`.
    Utility to convert old files to the new file as been added
    {meth}`~scvi.model.base.BaseModelClass.convert_legacy_save` ([#1505]).
- Breaking changes to {class}`~scvi.model.DestVI` as specified above ([#1457]).

#### Bug Fixes

- {meth}`~scvi.model.base.RNASeqMixin.get_likelihood_parameters` fix for `n_samples > 1` and
    `dispersion="gene_cell"` [#1504].
- Fix backwards compatibility for legacy TOTALVI models [#1502].

#### Contributors

- [@cane11]
- [@jjhong922]
- [@adamgayoso]
- [@romain-lopez]

## Version 0.15

### New in 0.15.5 (2022-04-06)

#### Changes

- Add common types file [#1467].
- New default is to not pin memory during training when using a GPU. This is much better for shared
    GPU environments without any performance regression [#1473].

#### Bug fixes

- Fix LDA user guide bugs [#1479].
- Fix unnecessary warnings, double logging [#1475].

#### Contributors

- [@jjhong922]
- [@adamgayoso]

### New in 0.15.4 (2022-03-28)

#### Changes

- Add peakVI publication reference [#1463].
- Update notebooks with new install functionality for Colab [#1466].
- Simplify changing the training plan for pyro [#1470].
- Optionally scale ELBO by a scalar in {class}`~scvi.train.PyroTrainingPlan` [#1469].

#### Bug fixes

#### Contributors

- [@jjhong922]
- [@adamgayoso]
- [@vitkl]

### New in 0.15.3 (2022-03-24)

#### Changes

#### Bug fixes

- Raise `NotImplementedError` when `categorical_covariate_keys` are used with
    {meth}`scvi.model.SCANVI.load_query_data`. ([#1458]).
- Fix behavior when `continuous_covariate_keys` are used with {meth}`scvi.model.SCANVI.classify`.
    ([#1458]).
- Unlabeled category values are automatically populated when
    {meth}`scvi.model.SCANVI.load_query_data` run on `adata_target` missing labels column.
    ([#1458]).
- Fix dataframe rendering in dark mode docs ([#1448])
- Fix variance constraint in {class}`~scvi.model.AmortizedLDA` that set an artifical bound on
    latent topic variance ([#1445]).
- Fix {meth}`scvi.model.base.ArchesMixin.prepare_query_data` to work cross device (e.g., model
    trained on cuda but method used on cpu; see [#1451]).

#### Contributors

- [@jjhong922]
- [@adamgayoso]

### New in 0.15.2 (2022-03-15)

#### Changes

- Remove setuptools pinned requirement due to new PyTorch 1.11 fix ([#1436]).
- Switch to myst-parsed markdown for docs ([#1435]).
- Add `prepare_query_data(adata, reference_model)` to {class}`~scvi.model.base.ArchesMixin` to
    enable query data cleaning prior to reference mapping ([#1441]).
- Add Human Lung Cell Atlas tutorial ([#1442]).

#### Bug fixes

- Errors when arbitrary kwargs are passed into `setup_anndata()` ([#1439]).
- Fix {class}`scvi.external.SOLO` to use `train_size=0.9` by default, which enables early stopping
    to work properly ([#1438]).
- Fix scArches version warning ([#1431]).
- Fix backwards compat for {class}`~scvi.model.SCANVI` loading ([#1441]).

#### Contributors

- [@jjhong922]
- [@adamgayoso]
- [@grst]

### New in 0.15.1 (2022-03-11)

#### Changes

- Remove `labels_key` from {class}`~scvi.model.MULTIVI` as it is not used in the model ([#1393]).
- Use scvi-tools mean/inv_disp parameterization of negative binomial for
    {class}`~scvi.model.JaxSCVI` likelihood ([#1386]).
- Use `setup` for Flax-based modules ([#1403]).
- Reimplement {class}`~scvi.module.JaxVAE` using inference/generative paradigm with
    {class}`~scvi.module.base.JaxBaseModuleClass` ([#1406]).
- Use multiple particles optionally in {class}`~scvi.model.JaxSCVI` ([#1385]).
- {class}`~scvi.external.SOLO` no longer warns about count data ([#1411]).
- Class docs are now one page on docs site ([#1415]).
- Copied AnnData objects are assigned a new uuid and transfer is attempted ([#1416]).

#### Bug fixes

- Fix an issue with using gene lists and proteins lists as well as `transform_batch` for
    {class}`~scvi.model.TOTALVI` ([#1413]).
- Error gracefully when NaNs present in {class}`~scvi.data.fields.CategoricalJointObsmField`
    ([#1417]).

#### Contributors

- [@jjhong922]
- [@adamgayoso]

### New in 0.15.0 (2022-02-28)

In this release, we have completely refactored the logic behind our data handling strategy (i.e.
`setup_anndata`) to allow for:

1. Readable data handling for existing models.
1. Modular code for easy addition of custom data fields to incorporate into models.
1. Avoidance of unexpected edge cases when more than one model is instantiated in one session.

**Important Note:** This change will not break pipelines for model users (with the exception of a
small change to {class}`~scvi.model.SCANVI`). However, there are several breaking changes for model
developers. The data handling tutorial goes over these changes in detail.

This refactor is centered around the new {class}`~scvi.data.AnnDataManager` class which
orchestrates any data processing necessary for scvi-tools and stores necessary information, rather
than adding additional fields to the AnnData input.

:::{figure} docs/\_static/img/anndata_manager_schematic.svg
:align: center
:alt: Schematic of data handling strategy with AnnDataManager
:class: img-fluid

Schematic of data handling strategy with {class}`~scvi.data.AnnDataManager`
:::

We also have an exciting new experimental Jax-based scVI implementation via
{class}`~scvi.model.JaxSCVI`. While this implementation has limited functionality, we have found it
to be substantially faster than the PyTorch-based implementation. For example, on a 10-core Intel
CPU, Jax on only a CPU can be as fast as PyTorch with a GPU (RTX3090). We will be planning further
Jax integrations in the next releases.

#### Changes

- Major refactor to data handling strategy with the introduction of
    {class}`~scvi.data.AnnDataManager` ([#1237]).
- Prevent clobbering between models using the same AnnData object with model instance specific
    {class}`~scvi.data.AnnDataManager` mappings ([#1342]).
- Add `size_factor_key` to {class}`~scvi.model.SCVI`, {class}`~scvi.model.MULTIVI`,
    {class}`~scvi.model.SCANVI`, and {class}`~scvi.model.TOTALVI` ([#1334]).
- Add references to the scvi-tools journal publication to the README ([#1338], [#1339]).
- Addition of {func}`scvi.model.utils.mde` ([#1372]) for accelerated visualization of scvi-tools
    embeddings.
- Documentation and user guide fixes ([#1364], [#1361])
- Fix for {class}`~scvi.external.SOLO` when {class}`~scvi.model.SCVI` was setup with a `labels_key`
    ([#1354])
- Updates to tutorials ([#1369], [#1371])
- Furo docs theme ([#1290])
- Add {class}`scvi.model.JaxSCVI` and {class}`scvi.module.JaxVAE`, drop Numba dependency for
    checking if data is count data ([#1367]).

#### Breaking changes

- The keyword argument `run_setup_anndata` has been removed from built-in datasets since there is
    no longer a model-agnostic `setup_anndata` method ([#1237]).

- The function `scvi.model._metrics.clustering_scores` has been removed due to incompatbility with
    new data handling ([#1237]).

- {class}`~scvi.model.SCANVI` now takes `unlabeled_category` as an argument to
    {meth}`~scvi.model.SCANVI.setup_anndata` rather than on initialization ([#1237]).

- `setup_anndata` is now a class method on model classes and requires specific function calls to
    ensure proper {class}`~scvi.data.AnnDataManager` setup and model save/load. Any model
    inheriting from {class}`~scvi.model.base.BaseModelClass` will need to re-implement this method
    ([#1237]).

    - To adapt existing custom models to v0.15.0, one can references the guidelines below. For
        some examples of how this was done for the existing models in the codebase, please
        reference the following PRs: ([#1301], [#1302]).
    - `scvi._CONSTANTS` has been changed to `scvi.REGISTRY_KEYS`.
    - `setup_anndata()` functions are now class functions and follow a specific structure. Please
        refer to {meth}`~scvi.model.SCVI.setup_anndata` for an example.
    - `scvi.data.get_from_registry()` has been removed. This method can be replaced by
        {meth}`scvi.data.AnnDataManager.get_from_registry`.
    - The setup dict stored directly on the AnnData object, `adata["_scvi"]`, has been deprecated.
        Instead, this information now lives in {attr}`scvi.data.AnnDataManager.registry`.
    - The data registry can be accessed at {attr}`scvi.data.AnnDataManager.data_registry`.
    - Summary stats can be accessed at {attr}`scvi.data.AnnDataManager.summary_stats`.
    - Any field-specific information (e.g. `adata.obs["categorical_mappings"]`) now lives in
        field-specific state registries. These can be retrieved via the function
        {meth}`~scvi.data.AnnDataManager.get_state_registry`.
    - `register_tensor_from_anndata()` has been removed. To register tensors with no relevant
        `AnnDataField` subclass, create a new a new subclass of
        {class}`~scvi.data.fields.BaseAnnDataField` and add it to appropriate model's
        `setup_anndata()` function.

#### Contributors

- [@jjhong922]
- [@adamgayoso]
- [@watiss]

## Version 0.14

### New in 0.14.6 (2021-02-05)

Bug fixes, minor improvements of docs, code formatting.

#### Changes

- Update black formatting to stable release ([#1324])
- Refresh readme, move tasks image to docs ([#1311]).
- Add 0.14.5 release note to index ([#1296]).
- Add test to ensure extra {class}`~scvi.model.SCANVI` training of a pre-trained
    {class}`~scvi.model.SCVI` model does not change original model weights ([#1284]).
- Fix issue in {class}`~scvi.model.TOTALVI` protein background prior initialization to not include
    protein measurements that are known to be missing ([#1282]).
- Upper bound setuptools due to PyTorch import bug ([#1309]).

#### Contributors

- [@adamgayoso]
- [@watiss]
- [@jjhong922]

### New in 0.14.5 (2021-11-22)

Bug fixes, new tutorials.

#### Changes

- Fix `kl_weight` floor for Pytorch-based models ([#1269]).
- Add support for more Pyro guides ([#1267]).
- Update scArches, harmonization tutorials, add basic R tutorial, tabula muris label transfer
    tutorial ([#1274]).

#### Contributors

- [@adamgayoso]
- [@jjhong922]
- [@watiss]
- [@vitkl]

### New in 0.14.4 (2021-11-16)

Bug fixes, some tutorial improvements.

#### Changes

- `kl_weight` handling for Pyro-based models ([#1242]).
- Allow override of missing protein inference in {class}`~scvi.model.TOTALVI` ([#1251]). This
    allows to treat all 0s in a particular batch for one protein as biologically valid.
- Fix load documentation (e.g., {meth}`~scvi.model.SCVI.load`, {meth}`~scvi.model.TOTALVI.load`)
    ([#1253]).
- Fix model history on load with Pyro-based models ([#1255]).
- Model construction tutorial uses new static setup anndata ([#1257]).
- Add codebase overview figure to docs ([#1231]).

#### Contributors

- [@adamgayoso]
- [@jjhong922]
- [@watiss]

### New in 0.14.3 (2021-10-19)

Bug fix.

#### Changes

- Bug fix to {func}`~scvi.model.base.BaseModelClass` to retain tensors registered by
    `register_tensor_from_anndata` ([#1235]).
- Expose an instance of our `DocstringProcessor` to aid in documenting derived implementations of
    `setup_anndata` method ([#1235]).

#### Contributors

- [@adamgayoso]
- [@jjhong922]
- [@watiss]

### New in 0.14.2 (2021-10-18)

Bug fix and new tutorial.

#### Changes

- Bug fix in {class}`~scvi.external.RNAStereoscope` where loss was computed with mean for a
    minibatch instead of sum. This ensures reproducibility with the original implementation ([#1228]).
- New Cell2location contributed tutorial ([#1232]).

#### Contributors

- [@adamgayoso]
- [@jjhong922]
- [@vitkl]
- [@watiss]

### New in 0.14.1 (2021-10-11)

Minor hotfixes.

#### Changes

- Filter out mitochrondrial genes as a preprocessing step in the Amortized LDA tutorial ([#1213])
- Remove `verbose=True` argument from early stopping callback ([#1216])

#### Contributors

- [@adamgayoso]
- [@jjhong922]
- [@watiss]

### New in 0.14.0 (2021-10-07)

In this release, we have completely revamped the scvi-tools documentation website by creating a
new set of user guides that provide:

1. The math behind each method (in a succinct, online methods-like way)
1. The relationship between the math and the functions associated with each model
1. The relationship between math variables and code variables

Our previous User Guide guide has been renamed to Tutorials and contains all of our existing
tutorials (including tutorials for developers).

Another noteworthy addition in this release is the implementation of the (amortized) Latent
Dirichlet Allocation (aka LDA) model applied to single-cell gene expression data. We have also
prepared a tutorial that demonstrates how to use this model, using a PBMC 10K dataset from 10x
Genomics as an example application.

Lastly, in this release we have made a change to reduce user and developer confusion by making the
previously global `setup_anndata` method a static class-specific method instead. This provides more
clarity on which parameters are applicable for this call, for each model class. Below is a
before/after for the DESTVI and TOTALVI model classes:

:::{figure} docs/\_static/img/setup_anndata_before_after.svg
:align: center
:alt: setup_anndata before and after
:class: img-fluid

`setup_anndata` before and after
:::

#### Changes

- Added fixes to support PyTorch Lightning 1.4 ([#1103])
- Simplified data handling in R tutorials with sceasy and addressed bugs in package installation
    ([#1122]).
- Moved library size distribution computation to model init ([#1123])
- Updated Contribution docs to describe how we backport patches ([#1129])
- Implemented Latent Dirichlet Allocation as a PyroModule ([#1132])
- Made `setup_anndata` a static method on model classes rather than one global function ([#1150])
- Used Pytorch Lightning's `seed_everything` method to set seed ([#1151])
- Fixed a bug in {class}`~scvi.model.base.PyroSampleMixin` for posterior sampling ([#1158])
- Added CITE-Seq datasets ([#1182])
- Added user guides to our documentation ([#1127], [#1157], [#1180], [#1193], [#1183], [#1204])
- Early stopping now prints the reason for stopping when applicable ([#1208])

#### Breaking changes

- `setup_anndata` is now an abstract method on model classes. Any model inheriting from
    {class}`~scvi.model.base.BaseModelClass` will need to implement this method ([#1150])

#### Contributors

- [@adamgayoso]
- [@PierreBoyeau]
- [@talashuach]
- [@jjhong922]
- [@watiss]
- [@mjayasur]
- [@vitkl]
- [@galenxing]

## Version 0.13

### New in 0.13.0 (2021-08-23)

#### Changes

- Added {class}`~scvi.model.MULTIVI` ([#1115], [#1118]).
- Documentation CSS tweaks ([#1116]).

#### Breaking changes

None!

#### Contributors

- [@adamgayoso]
- [@talashuach]
- [@jjhong922]

## Version 0.12

### New in 0.12.2 (2021-08-11)

#### Changes

- Updated `OrderedDict` typing import to support all Python 3.7 versions ([#1114]).

#### Breaking changes

None!

#### Contributors

- [@adamgayoso]
- [@galenxing]
- [@jjhong922]

### New in 0.12.1 (2021-07-29)

#### Changes

- Update Pytorch Lightning version dependency to `>=1.3,<1.4` ([#1104]).

#### Breaking changes

None!

#### Contributors

- [@adamgayoso]
- [@galenxing]

### New in 0.12.0 (2021-07-15)

This release adds features for tighter integration with Pyro for model development, fixes for
{class}`~scvi.external.SOLO`, and other enhancements. Users of {class}`~scvi.external.SOLO` are
strongly encouraged to upgrade as previous bugs will affect performance.

#### Enchancements

- Add {class}`scvi.model.base.PyroSampleMixin` for easier posterior sampling with Pyro ([#1059]).
- Add {class}`scvi.model.base.PyroSviTrainMixin` for automated training of Pyro models ([#1059]).
- Ability to pass kwargs to {class}`~scvi.module.Classifier` when using
    {class}`~scvi.external.SOLO` ([#1078]).
- Ability to get doublet predictions for simulated doublets in {class}`~scvi.external.SOLO`
    ([#1076]).
- Add "comparison" column to differential expression results ([#1074]).
- Clarify {class}`~scvi.external.CellAssign` size factor usage. See class docstring.

#### Changes

- Update minimum Python version to `3.7.2` ([#1082]).
- Slight interface changes to {class}`~scvi.train.PyroTrainingPlan`. `"elbo_train"` and
    `"elbo_test"` are now the average over minibatches as ELBO should be on scale of full data and
    `optim_kwargs` can be set on initialization of training plan ([#1059], [#1101]).
- Use pandas read pickle function for pbmc dataset metadata loading ([#1099]).
- Adds `n_samples_overall` parameter to functions for denoised expression/accesibility/etc. This is
    used in during differential expression ([#1090]).
- Ignore configure optimizers warning when training Pyro-based models ([#1064]).

#### Bug fixes

- Fix scale of library size for simulated doublets and expression in {class}`~scvi.external.SOLO`
    when using observed library size to train original {class}`~scvi.model.SCVI` model ([#1078],
    [#1085]). Currently, library sizes in this case are not appropriately put on the log scale.
- Fix issue where anndata setup with a layer led to errors in {class}`~scvi.external.SOLO`
    ([#1098]).
- Fix `adata` parameter of {func}`scvi.external.SOLO.from_scvi_model`, which previously did nothing
    ([#1078]).
- Fix default `max_epochs` of {class}`~scvi.model.SCANVI` when initializing using pre-trained model
    of {class}`~scvi.model.SCVI` ([#1079]).
- Fix bug in `predict()` function of {class}`~scvi.model.SCANVI`, which only occurred for soft
    predictions ([#1100]).

#### Breaking changes

None!

#### Contributors

- [@vitkl]
- [@adamgayoso]
- [@galenxing]
- [@PierreBoyeau]
- [@Munfred]
- [@njbernstein]
- [@mjayasur]

## Version 0.11

### New in 0.11.0 (2021-05-23)

From the user perspective, this release features the new differential expression functionality (to
be described in a manuscript). For now, it is accessible from
{func}`~scvi.model.SCVI.differential_expression`. From the developer perspective, we made changes
with respect to {class}`scvi.dataloaders.DataSplitter` and surrounding the Pyro backend. Finally,
we also made changes to adapt our code to PyTorch Lightning version 1.3.

#### Changes

- Pass `n_labels` to {class}`~scvi.module.VAE` from {class}`~scvi.model.SCVI` ([#1055]).
- Require PyTorch lightning > 1.3, add relevant fixes ([#1054]).
- Add DestVI reference ([#1060]).
- Add PeakVI links to README ([#1046]).
- Automatic delta and eps computation in differential expression ([#1043]).
- Allow doublet ratio parameter to be changed for used in SOLO ([#1066]).

#### Bug fixes

- Fix an issue where `transform_batch` options in {class}`~scvi.model.TOTALVI` was accidentally
    altering the batch encoding in the encoder, which leads to poor results ([#1072]). This bug was
    introduced in version 0.9.0.

#### Breaking changes

These breaking changes do not affect the user API; though will impact model developers.

- Use PyTorch Lightning data modules for {class}`scvi.dataloaders.DataSplitter` ([#1061]). This
    induces a breaking change in the way the data splitter is used. It is no longer callable and
    now has a `setup` method. See {class}`~scvi.train.TrainRunner` and its source code, which is
    straightforward.
- No longer require training plans to be initialized with `n_obs_training` argument ([#1061]).
    `n_obs_training` is now a property that can be set before actual training to rescale the loss.
- Log Pyro loss as `train_elbo` and sum over steps ([#1071])

#### Contributors

- [@adamgayoso]
- [@romain-lopez]
- [@PierreBoyeau]
- [@talashuach]
- [@cataclysmus]
- [@njbernstein]

## Version 0.10

### New in 0.10.1 (2021-05-04)

#### Changes

- Includes new optional variance parameterization for the `Encoder` module ([#1037]).
- Provides new way to select subpopulations for DE using Pandas queries ([#1041]).
- Update reference to peakVI ([#1046]).
- Pin Pytorch Lightning version to \<1.3

#### Contributors

- [@adamgayoso]
- [@PierreBoyeau]
- [@talashuach]

### New in 0.10.0 (2021-04-20)

#### Changes

- PeakVI minor enhancements to differential accessibility and fix scArches support ([#1019])
- Add DestVI to the codebase ([#1011])
- Versioned tutorial links ([#1005])
- Remove old VAEC ([#1006])
- Use `.numpy()` to convert torch tensors to numpy ndarrays ([#1016])
- Support backed AnnData ([#1017]), just load anndata with `scvi.data.read_h5ad(path, backed='r+')`
- Solo interface enhancements ([#1009])
- Updated README ([#1028])
- Use Python warnings instead of logger warnings ([#1021])
- Change totalVI protein background default to `False` is fewer than 10 proteins used ([#1034])

#### Bug fixes

- Fix `SaveBestState` warning ([#1024])
- New default SCANVI max epochs if loaded with pretrained SCVI model ([#1025]), restores old
    `<v0.9` behavior.
- Fix marginal log likelihood computation, which was only being computed on final minibatch of a
    dataloader. This bug was introduced in the `0.9.X` versions ([#1033]).
- Fix bug where extra categoricals were not properly extended in `transfer_anndata_setup` ([#1030]).

#### Contributors

- [@adamgayoso]
- [@romain-lopez]
- [@talashuach]
- [@mjayasur]
- [@wukathy]
- [@PierreBoyeau]
- [@morris-frank]

## Version 0.9

### New in 0.9.1 (2021-03-20)

#### Changes

- Update Pyro module backend to better enfore usage of `model` and `guide`, automate passing of
    number of training examples to Pyro modules ([#990])
- Minimum Pyro version bumped ([#988])
- Improve docs clarity ([#989])
- Add glossary to developer user guide ([#999])
- Add num threads config option to `scvi.settings` ([#1001])
- Add CellAssign tutorial ([#1004])

#### Contributors

- [@adamgayoso]
- [@galenxing]
- [@mjayasur]
- [@wukathy]

### New in 0.9.0 (2021-03-03)

This release features our new software development kit for building new probabilistic models. Our
hope is that others will be able to develop new models by importing scvi-tools into their own
packages.

#### Important changes

From the user perspective, there are two package-wide API breaking changes and one
{class}`~scvi.model.SCANVI` specific breaking change enumerated below. From the method developer
perspective, the entire model backend has been revamped using PyTorch Lightning, and no old code
will be compatible with this and future versions. Also, we dropped support for Python 3.6.

##### Breaking change: The `train` method

- `n_epochs` is now `max_epochs` for consistency with PytorchLightning and to better relect the
    functionality of the parameter.
- `use_cuda` is now `use_gpu` for consistency with PytorchLightning.
- `frequency` is now `check_val_every_n_epoch` for consistency with PytorchLightning.
- `train_fun_kwargs` and `kwargs` throughout the `train()` methods in the codebase have been
    removed and various arguments have been reorganized into `plan_kwargs` and `trainer_kwargs`.
    Generally speaking, `plan_kwargs` deal with model optimization like kl warmup, while
    `trainer_kwargs` deal with the actual training loop like early stopping.

##### Breaking change: GPU handling

- `use_cuda` was removed from the init of each model and was not replaced by `use_gpu`. By default
    every model is intialized on CPU but can be moved to a device via `model.to_device()`. If a
    model is trained with `use_gpu=True` the model will remain on the GPU after training.
- When loading saved models, scvi-tools will always attempt to load the model on GPU unless
    otherwise specified.
- We now support specifying which GPU device to use if there are multiple available GPUs.

##### Breaking change: {class}`~scvi.model.SCANVI`

- {class}`~scvi.model.SCANVI` no longer pretrains an {class}`~scvi.model.SCVI` model by default.
    This functionality however is preserved via the new {func}`~scvi.model.SCANVI.from_scvi_model`
    method.
- `n_epochs_unsupervised` and `n_epochs_semisupervised` have been removed from `train`. It has been
    replaced with `max_epochs` for semisupervised training.
- `n_samples_per_label` is a new argument which will subsample the number of labelled training
    examples to train on per label each epoch.

#### New Model Implementations

- {class}`~scvi.model.PEAKVI` implementation ([#877], [#921])
- {class}`~scvi.external.SOLO` implementation ([#923], [#933])
- {class}`~scvi.external.CellAssign` implementation ([#940])
- {class}`~scvi.external.RNAStereoscope` and {class}`~scvi.external.SpatialStereoscope`
    implementation ([#889], [#959])
- Pyro integration via {class}`~scvi.module.base.PyroBaseModuleClass` ([#895] [#903], [#927],
    [#931])

#### Enhancements

- {class}`~scvi.model.SCANVI` bug fixes ([#879])
- {class}`~scvi.external.GIMVI` moved to external api ([#885])
- {class}`~scvi.model.TOTALVI`, {class}`~scvi.model.SCVI`, and {class}`~scvi.model.SCANVI` now
    support multiple covariates ([#886])
- Added callback for saving the best state of a model ([#887])
- Option to disable progress bar ([#905])
- load() documentation improvements ([#913])
- updated tutorials, guides, documentation ([#924], [#925], [#929], [#934], [#947], [#971])
- track is now public ([#938])
- {class}`~scvi.model.SCANVI` now logs classficiation loss ([#966])
- get_likelihood_parameter() bug ([#967])
- model.history are now pandas DataFrames ([#949])

#### Contributors

- [@adamgayoso]
- [@galenxing]
- [@romain-lopez]
- [@wukathy]
- [@giovp]
- [@njbernstein]
- [@saketkc]

## Version 0.8

### New in 0.8.1 (2020-12-23)

#### Enhancements

- `freeze_classifier` option in {func}`~scvi.model.SCANVI.load_query_data` for the case when
- `weight_decay` passed to {func}`~scvi.model.SCANVI.train` also passes to `ClassifierTrainer`

### New in 0.8.0 (2020-12-17)

#### Enhancements

##### Online updates of {class}`~scvi.model.SCVI`, {class}`~scvi.model.SCANVI`, and {class}`~scvi.model.TOTALVI` with the scArches method <!-- markdownlint-disable -->

It is now possible to iteratively update these models with new samples, without altering the model
for the "reference" population. Here we use the
[scArches method](https://github.com/theislab/scarches). For usage, please see the tutorial in the
user guide.

To enable scArches in our models, we added a few new options. The first is `encode_covariates`,
which is an `SCVI` option to encode the one-hotted batch covariate. We also allow users to exchange
batch norm in the encoder and decoder with layer norm, which can be though of as batch norm but per
cell. As the layer norm we use has no parameters, it's a bit faster than models with batch norm. We
don't find many differences between using batch norm or layer norm in our models, though we have
kept defaults the same in this case. To run scArches effectively, batch norm should be exhanged
with layer norm.

##### Empirical initialization of protein background parameters with totalVI

The learned prior parameters for the protein background were randomly initialized. Now, they can be
set with the `empirical_protein_background_prior` option in {class}`~scvi.model.TOTALVI`. This
option fits a two-component Gaussian mixture model per cell, separating those proteins that are
background for the cell and those that are foreground, and aggregates the learned mean and variance
of the smaller component across cells. This computation is done per batch, if the `batch_key` was
registered. We emphasize this is just for the initialization of a learned parameter in the model.

##### Use observed library size option

Many of our models like `SCVI`, `SCANVI`, and {class}`~scvi.model.TOTALVI` learn a latent library
size variable. The option `use_observed_lib_size` may now be passed on model initialization. We
have set this as `True` by default, as we see no regression in performance, and training is a bit
faster.

#### Important changes

- To facilitate these enhancements, saved {class}`~scvi.model.TOTALVI` models from previous
    versions will not load properly. This is due to an architecture change of the totalVI encoder,
    related to latent library size handling.
- The default latent distribtuion for {class}`~scvi.model.TOTALVI` is now `"normal"`.
- Autotune was removed from this release. We could not maintain the code given the new API changes
    and we will soon have alternative ways to tune hyperparameters.
- Protein names during `setup_anndata` are now stored in `adata.uns["_scvi"]["protein_names"]`,
    instead of `adata.uns["scvi_protein_names"]`.

#### Bug fixes

- Fixed an issue where the unlabeled category affected the SCANVI architecture prior distribution.
    Unfortunately, by fixing this bug, loading previously trained (\<v0.8.0)
    {class}`~scvi.model.SCANVI` models will fail.

## Version 0.7

### New in 0.7.1 (2020-10-20)

This small update provides access to our new Discourse forum from the documentation.

### New in 0.7.0 (2020-10-14)

scvi is now scvi-tools. Version 0.7 introduces many breaking changes. The best way to learn how to
use scvi-tools is with our documentation and tutorials.

- New high-level API and data loading, please see tutorials and examples for usage.
- `GeneExpressionDataset` and associated classes have been removed.
- Built-in datasets now return `AnnData` objects.
- `scvi-tools` now relies entirely on the [AnnData] format.
- `scvi.models` has been moved to `scvi.core.module`.
- `Posterior` classes have been reduced to wrappers on `DataLoaders`
- `scvi.inference` has been split to `scvi.core.data_loaders` for `AnnDataLoader` classes and
    `scvi.core.trainers` for trainer classes.
- Usage of classes like `Trainer` and `AnnDataLoader` now require the `AnnData` data object as
    input.

## Pre-Version 0.7

### scvi History

The scvi-tools package used to be scvi. This page commemorates all the hard work on the scvi
package by our numerous contributors.

#### Contributors

- [@romain]
- [@adam]
- [@eddie]
- [@jeff]
- [@pierre]
- [@max]
- [@yining]
- [@gabriel]
- [@achille]
- [@chenling]
- [@jules]
- [@david-kelley]
- [@william-yang]
- [@oscar]
- [@casey-greene]
- [@jamie-morton]
- [@valentine-svensson]
- [@stephen-flemming]
- [@michael-raevsky]
- [@james-webber]
- [@galen]
- [@francesco-brundu]
- [@primoz-godec]
- [@eduardo-beltrame]
- [@john-reid]
- [@han-yuan]
- [@gokcen-eraslan]

#### 0.6.7 (2020-8-05)

- downgrade anndata>=0.7 and scanpy>=1.4.6 [@galen]
- make loompy optional, raise sckmisc import error [@adam]
- fix PBMCDataset download bug [@galen]
- fix AnnDatasetFromAnnData \_X in adata.obs bug [@galen]

#### 0.6.6 (2020-7-08)

- add tqdm to within cluster DE genes [@adam]
- restore tqdm to use simple bar instead of ipywidget [@adam]
- move to numpydoc for doctstrings [@adam]
- update issues templates [@adam]
- Poisson variable gene selection [@valentine-svensson]
- BrainSmallDataset set defualt save_path_10X [@gokcen-eraslan]
- train_size must be float between 0.0 and 1.0 [@galen]
- bump dependency versions [@galen]
- remove reproducibility notebook [@galen]
- fix scanVI dataloading [@pierre]

#### 0.6.5 (2020-5-10)

- updates to totalVI posterior functions and notebooks [@adam]
- update seurat v3 HVG selection now using skmisc loess [@adam]

#### 0.6.4 (2020-4-14)

- add back Python 3.6 support [@adam]
- get_sample_scale() allows gene selection [@valentine-svensson]
- bug fix to the dataset to anndata method with how cell measurements are stored [@adam]
- fix requirements [@adam]

#### 0.6.3 (2020-4-01)

- bug in version for Louvian in setup.py [@adam]

#### 0.6.2 (2020-4-01)

- update highly variable gene selection to handle sparse matrices [@adam]
- update DE docstrings [@pierre]
- improve posterior save load to also handle subclasses [@pierre]
- Create NB and ZINB distributions with torch and refactor code accordingly [@pierre]
- typos in autozivae [@achille]
- bug in csc sparse matrices in anndata data loader [@adam]

#### 0.6.1 (2020-3-13)

- handles gene and cell attributes with the same name [@han-yuan]
- fixes anndata overwriting when loading [@adam], [@pierre]
- formatting in basic tutorial [@adam]

#### 0.6.0 (2020-2-28)

- updates on TotalVI and LDVAE [@adam]
- fix documentation, compatibility and diverse bugs [@adam], [@pierre] [@romain]
- fix for external module on scanpy [@galen]

#### 0.5.0 (2019-10-17)

- do not automatically upper case genes [@adam]
- AutoZI [@oscar]
- Made the intro tutorial more user friendly [@adam]
- Tests for LDVAE notebook [@adam]
- black codebase [@achille] [@gabriel] [@adam]
- fix compatibility issues with sklearn and numba [@romain]
- fix Anndata [@francesco-brundu]
- docstring, totalVI, totalVI notebook and CITE-seq data [@adam]
- fix type [@eduardo-beltrame]
- fixing installation guide [@jeff]
- improved error message for dispersion [@stephen-flemming]

#### 0.4.1 (2019-08-03)

- docstring [@achille]
- differential expression [@oscar] [@pierre]

#### 0.4.0 (2019-07-25)

- gimVI [@achille]
- synthetic correlated datasets, fixed bug in marginal log likelihood [@oscar]
- autotune, dataset enhancements [@gabriel]
- documentation [@jeff]
- more consistent posterior API, docstring, validation set [@adam]
- fix anndataset [@michael-raevsky]
- linearly decoded VAE [@valentine-svensson]
- support for scanpy, fixed bugs, dataset enhancements [@achille]
- fix filtering bug, synthetic correlated datasets, docstring, differential expression [@pierre]
- better docstring [@jamie-morton]
- classifier based on library size for doublet detection [@david-kelley]

#### 0.3.0 (2019-05-03)

- corrected notebook [@jules]
- added UMAP and updated harmonization code [@chenling] [@romain]
- support for batch indices in csvdataset [@primoz-godec]
- speeding up likelihood computations [@william-yang]
- better anndata interop [@casey-greene]
- early stopping based on classifier accuracy [@david-kelley]

#### 0.2.4 (2018-12-20)

- updated to torch v1 [@jules]
- added stress tests for harmonization [@chenling]
- fixed autograd breaking [@romain]
- make removal of empty cells more efficient [@john-reid]
- switch to os.path.join [@casey-greene]

#### 0.2.2 (2018-11-08)

- added baselines and datasets for sMFISH imputation [@jules]
- added harmonization content [@chenling]
- fixing bugs on DE [@romain]

#### 0.2.0 (2018-09-04)

- annotation notebook [@eddie]
- Memory footprint management [@jeff]
- updated early stopping [@max]
- docstring [@james-webber]

#### 0.1.6 (2018-08-08)

- MMD and adversarial inference wrapper [@eddie]
- Documentation [@jeff]
- smFISH data imputation [@max]

#### 0.1.5 (2018-07-24)

- Dataset additions [@eddie]
- Documentation [@yining]
- updated early stopping [@max]

#### 0.1.3 (2018-06-22)

- Notebook enhancement [@yining]
- Semi-supervision [@eddie]

#### 0.1.2 (2018-06-13)

- First release on PyPi
- Skeleton code & dependencies [@jeff]
- Unit tests [@max]
- PyTorch implementation of scVI [@eddie] [@max]
- Dataset preprocessing [@eddie] [@max] [@yining]

#### 0.1.0 (2017-09-05)

- First scVI TensorFlow version [@romain]

[#1001]: https://github.com/YosefLab/scvi-tools/pull/1001
[#1004]: https://github.com/YosefLab/scvi-tools/pull/1004
[#1005]: https://github.com/YosefLab/scvi-tools/pull/1005
[#1006]: https://github.com/YosefLab/scvi-tools/pull/1006
[#1009]: https://github.com/YosefLab/scvi-tools/pull/1009
[#1011]: https://github.com/YosefLab/scvi-tools/pull/1011
[#1016]: https://github.com/YosefLab/scvi-tools/pull/1016
[#1017]: https://github.com/YosefLab/scvi-tools/pull/1017
[#1019]: https://github.com/YosefLab/scvi-tools/pull/1019
[#1021]: https://github.com/YosefLab/scvi-tools/pull/1021
[#1024]: https://github.com/YosefLab/scvi-tools/pull/1025
[#1025]: https://github.com/YosefLab/scvi-tools/pull/1025
[#1028]: https://github.com/YosefLab/scvi-tools/pull/1028
[#1030]: https://github.com/YosefLab/scvi-tools/pull/1033
[#1033]: https://github.com/YosefLab/scvi-tools/pull/1033
[#1034]: https://github.com/YosefLab/scvi-tools/pull/1034
[#1037]: https://github.com/YosefLab/scvi-tools/pull/1037
[#1041]: https://github.com/YosefLab/scvi-tools/pull/1041
[#1043]: https://github.com/YosefLab/scvi-tools/pull/1043
[#1046]: https://github.com/YosefLab/scvi-tools/pull/1046
[#1054]: https://github.com/YosefLab/scvi-tools/pull/1054
[#1055]: https://github.com/YosefLab/scvi-tools/pull/1055
[#1059]: https://github.com/YosefLab/scvi-tools/pull/1059
[#1060]: https://github.com/YosefLab/scvi-tools/pull/1060
[#1061]: https://github.com/YosefLab/scvi-tools/pull/1061
[#1064]: https://github.com/YosefLab/scvi-tools/pull/1064
[#1066]: https://github.com/YosefLab/scvi-tools/pull/1066
[#1071]: https://github.com/YosefLab/scvi-tools/pull/1071
[#1072]: https://github.com/YosefLab/scvi-tools/pull/1072
[#1074]: https://github.com/YosefLab/scvi-tools/pull/1074
[#1076]: https://github.com/YosefLab/scvi-tools/pull/1076
[#1078]: https://github.com/YosefLab/scvi-tools/pull/1078
[#1079]: https://github.com/YosefLab/scvi-tools/pull/1079
[#1082]: https://github.com/YosefLab/scvi-tools/pull/1082
[#1085]: https://github.com/YosefLab/scvi-tools/pull/1085
[#1090]: https://github.com/YosefLab/scvi-tools/pull/1090
[#1098]: https://github.com/YosefLab/scvi-tools/pull/1098
[#1099]: https://github.com/YosefLab/scvi-tools/pull/1099
[#1100]: https://github.com/YosefLab/scvi-tools/pull/1100
[#1101]: https://github.com/YosefLab/scvi-tools/pull/1101
[#1103]: https://github.com/YosefLab/scvi-tools/pull/1103
[#1104]: https://github.com/YosefLab/scvi-tools/pull/1104
[#1114]: https://github.com/YosefLab/scvi-tools/pull/1114
[#1115]: https://github.com/YosefLab/scvi-tools/pull/1115
[#1116]: https://github.com/YosefLab/scvi-tools/pull/1116
[#1118]: https://github.com/YosefLab/scvi-tools/pull/1118
[#1122]: https://github.com/YosefLab/scvi-tools/pull/1122
[#1123]: https://github.com/YosefLab/scvi-tools/pull/1123
[#1127]: https://github.com/YosefLab/scvi-tools/pull/1127
[#1129]: https://github.com/YosefLab/scvi-tools/pull/1129
[#1132]: https://github.com/YosefLab/scvi-tools/pull/1132
[#1150]: https://github.com/YosefLab/scvi-tools/pull/1150
[#1151]: https://github.com/YosefLab/scvi-tools/pull/1151
[#1157]: https://github.com/YosefLab/scvi-tools/pull/1157
[#1158]: https://github.com/YosefLab/scvi-tools/pull/1158
[#1180]: https://github.com/YosefLab/scvi-tools/pull/1180
[#1182]: https://github.com/YosefLab/scvi-tools/pull/1182
[#1183]: https://github.com/YosefLab/scvi-tools/pull/1183
[#1193]: https://github.com/YosefLab/scvi-tools/pull/1193
[#1204]: https://github.com/YosefLab/scvi-tools/pull/1204
[#1208]: https://github.com/YosefLab/scvi-tools/pull/1208
[#1213]: https://github.com/YosefLab/scvi-tools/pull/1213
[#1216]: https://github.com/YosefLab/scvi-tools/pull/1216
[#1228]: https://github.com/YosefLab/scvi-tools/pull/1228
[#1231]: https://github.com/YosefLab/scvi-tools/pull/1231
[#1232]: https://github.com/YosefLab/scvi-tools/pull/1232
[#1235]: https://github.com/YosefLab/scvi-tools/pull/1235
[#1237]: https://github.com/YosefLab/scvi-tools/pull/1237
[#1242]: https://github.com/YosefLab/scvi-tools/pull/1242
[#1251]: https://github.com/YosefLab/scvi-tools/pull/1251
[#1253]: https://github.com/YosefLab/scvi-tools/pull/1253
[#1255]: https://github.com/YosefLab/scvi-tools/pull/1255
[#1257]: https://github.com/YosefLab/scvi-tools/pull/1257
[#1267]: https://github.com/YosefLab/scvi-tools/pull/1267
[#1269]: https://github.com/YosefLab/scvi-tools/pull/1269
[#1274]: https://github.com/YosefLab/scvi-tools/pull/1274
[#1282]: https://github.com/YosefLab/scvi-tools/pull/1282
[#1284]: https://github.com/YosefLab/scvi-tools/pull/1284
[#1290]: https://github.com/YosefLab/scvi-tools/pull/1290
[#1296]: https://github.com/YosefLab/scvi-tools/pull/1296
[#1301]: https://github.com/YosefLab/scvi-tools/pull/1301
[#1302]: https://github.com/YosefLab/scvi-tools/pull/1302
[#1309]: https://github.com/YosefLab/scvi-tools/pull/1309
[#1311]: https://github.com/YosefLab/scvi-tools/pull/1311
[#1324]: https://github.com/YosefLab/scvi-tools/pull/1324
[#1334]: https://github.com/YosefLab/scvi-tools/pull/1334
[#1338]: https://github.com/YosefLab/scvi-tools/pull/1338
[#1339]: https://github.com/YosefLab/scvi-tools/pull/1339
[#1342]: https://github.com/YosefLab/scvi-tools/pull/1342
[#1354]: https://github.com/YosefLab/scvi-tools/pull/1354
[#1356]: https://github.com/YosefLab/scvi-tools/pull/1356
[#1361]: https://github.com/YosefLab/scvi-tools/pull/1361
[#1364]: https://github.com/YosefLab/scvi-tools/pull/1364
[#1367]: https://github.com/YosefLab/scvi-tools/pull/1367
[#1369]: https://github.com/YosefLab/scvi-tools/pull/1369
[#1371]: https://github.com/YosefLab/scvi-tools/pull/1371
[#1372]: https://github.com/YosefLab/scvi-tools/pull/1372
[#1385]: https://github.com/YosefLab/scvi-tools/pull/1385
[#1386]: https://github.com/YosefLab/scvi-tools/pull/1386
[#1393]: https://github.com/YosefLab/scvi-tools/pull/1393
[#1403]: https://github.com/YosefLab/scvi-tools/pull/1403
[#1406]: https://github.com/YosefLab/scvi-tools/pull/1406
[#1411]: https://github.com/YosefLab/scvi-tools/pull/1411
[#1413]: https://github.com/YosefLab/scvi-tools/pull/1413
[#1415]: https://github.com/YosefLab/scvi-tools/pull/1415
[#1416]: https://github.com/YosefLab/scvi-tools/pull/1416
[#1417]: https://github.com/YosefLab/scvi-tools/pull/1417
[#1431]: https://github.com/YosefLab/scvi-tools/pull/1431
[#1435]: https://github.com/YosefLab/scvi-tools/pull/1435
[#1436]: https://github.com/YosefLab/scvi-tools/pull/1436
[#1438]: https://github.com/YosefLab/scvi-tools/pull/1438
[#1439]: https://github.com/YosefLab/scvi-tools/pull/1439
[#1441]: https://github.com/YosefLab/scvi-tools/pull/1441
[#1442]: https://github.com/YosefLab/scvi-tools/pull/1442
[#1445]: https://github.com/YosefLab/scvi-tools/pull/1445
[#1448]: https://github.com/YosefLab/scvi-tools/pull/1448
[#1451]: https://github.com/YosefLab/scvi-tools/pull/1451
[#1457]: https://github.com/YosefLab/scvi-tools/pull/1457
[#1458]: https://github.com/YosefLab/scvi-tools/pull/1458
[#1463]: https://github.com/YosefLab/scvi-tools/pull/1463
[#1466]: https://github.com/YosefLab/scvi-tools/pull/1466
[#1467]: https://github.com/YosefLab/scvi-tools/pull/1467
[#1469]: https://github.com/YosefLab/scvi-tools/pull/1469
[#1470]: https://github.com/YosefLab/scvi-tools/pull/1470
[#1473]: https://github.com/YosefLab/scvi-tools/pull/1473
[#1474]: https://github.com/YosefLab/scvi-tools/pull/1474
[#1475]: https://github.com/YosefLab/scvi-tools/pull/1475
[#1479]: https://github.com/YosefLab/scvi-tools/pull/1479
[#1498]: https://github.com/YosefLab/scvi-tools/pull/1498
[#1499]: https://github.com/YosefLab/scvi-tools/pull/1499
[#1501]: https://github.com/YosefLab/scvi-tools/pull/1501
[#1502]: https://github.com/YosefLab/scvi-tools/pull/1502
[#1504]: https://github.com/YosefLab/scvi-tools/pull/1504
[#1505]: https://github.com/YosefLab/scvi-tools/pull/1505
[#1506]: https://github.com/YosefLab/scvi-tools/pull/1506
[#1508]: https://github.com/YosefLab/scvi-tools/pull/1508
[#1515]: https://github.com/YosefLab/scvi-tools/pull/1515
[#1519]: https://github.com/YosefLab/scvi-tools/pull/1519
[#1520]: https://github.com/YosefLab/scvi-tools/pull/1520
[#1527]: https://github.com/YosefLab/scvi-tools/pull/1527
[#1529]: https://github.com/YosefLab/scvi-tools/pull/1529
[#1532]: https://github.com/YosefLab/scvi-tools/pull/1532
[#1542]: https://github.com/YosefLab/scvi-tools/pull/1542
[#1548]: https://github.com/YosefLab/scvi-tools/pull/1548
[#1555]: https://github.com/YosefLab/scvi-tools/pull/1555
[#1556]: https://github.com/YosefLab/scvi-tools/pull/1556
[#1566]: https://github.com/scverse/scvi-tools/issues/1566
[#1575]: https://github.com/YosefLab/scvi-tools/pull/1575
[#1580]: https://github.com/scverse/scvi-tools/pull/1580
[#1585]: https://github.com/YosefLab/scvi-tools/pull/1585
[#1595]: https://github.com/scverse/scvi-tools/pull/1595
[#1617]: https://github.com/scverse/scvi-tools/pull/1617
[#1618]: https://github.com/scverse/scvi-tools/pull/1618
[#1622]: https://github.com/scverse/scvi-tools/pull/1622
[#1629]: https://github.com/scverse/scvi-tools/pull/1629
[#1637]: https://github.com/scverse/scvi-tools/pull/1637
[#1639]: https://github.com/scverse/scvi-tools/pull/1639
[#1645]: https://github.com/scverse/scvi-tools/pull/1645
[#1657]: https://github.com/scverse/scvi-tools/pull/1657
[#1660]: https://github.com/scverse/scvi-tools/pull/1660
[#1665]: https://github.com/scverse/scvi-tools/pull/1665
[#1667]: https://github.com/scverse/scvi-tools/pull/1667
[#1671]: https://github.com/scverse/scvi-tools/pull/1671
[#1672]: https://github.com/YosefLab/scvi-tools/pull/1672
[#1674]: https://github.com/scverse/scvi-tools/pull/1674
[#1678]: https://github.com/scverse/scvi-tools/pull/1678
[#1683]: https://github.com/YosefLab/scvi-tools/pull/1683
[#1686]: https://github.com/scverse/scvi-tools/pull/1686
[#1689]: https://github.com/YosefLab/scvi-tools/pull/1689
[#1692]: https://github.com/scverse/scvi-tools/pull/1692
[#1695]: https://github.com/YosefLab/scvi-tools/pull/1695
[#1696]: https://github.com/YosefLab/scvi-tools/pull/1696
[#1697]: https://github.com/YosefLab/scvi-tools/pull/1697
[#1700]: https://github.com/scverse/scvi-tools/pull/1700
[#1702]: https://github.com/scverse/scvi-tools/pull/1702
[#1709]: https://github.com/YosefLab/scvi-tools/pull/1709
[#1710]: https://github.com/YosefLab/scvi-tools/pull/1710
[#1711]: https://github.com/YosefLab/scvi-tools/pull/1711
[#1719]: https://github.com/YosefLab/scvi-tools/pull/1719
[#1731]: https://github.com/YosefLab/scvi-tools/pull/1731
[#1732]: https://github.com/YosefLab/scvi-tools/pull/1732
[#1733]: https://github.com/YosefLab/scvi-tools/pull/1733
[#1737]: https://github.com/YosefLab/scvi-tools/pull/1737
[#1741]: https://github.com/YosefLab/scvi-tools/pull/1741
[#1743]: https://github.com/YosefLab/scvi-tools/pull/1743
[#1747]: https://github.com/YosefLab/scvi-tools/pull/1747
[#1749]: https://github.com/YosefLab/scvi-tools/pull/1749
[#1751]: https://github.com/YosefLab/scvi-tools/pull/1751
[#1773]: https://github.com/YosefLab/scvi-tools/pull/1773
[#877]: https://github.com/YosefLab/scvi-tools/pull/887
[#879]: https://github.com/YosefLab/scvi-tools/pull/879
[#885]: https://github.com/YosefLab/scvi-tools/pull/885
[#886]: https://github.com/YosefLab/scvi-tools/pull/886
[#887]: https://github.com/YosefLab/scvi-tools/pull/887
[#889]: https://github.com/YosefLab/scvi-tools/pull/889
[#895]: https://github.com/YosefLab/scvi-tools/pull/895
[#903]: https://github.com/YosefLab/scvi-tools/pull/903
[#905]: https://github.com/YosefLab/scvi-tools/pull/905
[#913]: https://github.com/YosefLab/scvi-tools/pull/913
[#921]: https://github.com/YosefLab/scvi-tools/pull/921
[#923]: https://github.com/YosefLab/scvi-tools/pull/923
[#924]: https://github.com/YosefLab/scvi-tools/pull/924
[#925]: https://github.com/YosefLab/scvi-tools/pull/925
[#927]: https://github.com/YosefLab/scvi-tools/pull/927
[#929]: https://github.com/YosefLab/scvi-tools/pull/929
[#931]: https://github.com/YosefLab/scvi-tools/pull/931
[#933]: https://github.com/YosefLab/scvi-tools/pull/933
[#934]: https://github.com/YosefLab/scvi-tools/pull/934
[#938]: https://github.com/YosefLab/scvi-tools/pull/938
[#940]: https://github.com/YosefLab/scvi-tools/pull/940
[#947]: https://github.com/YosefLab/scvi-tools/pull/947
[#949]: https://github.com/YosefLab/scvi-tools/pull/949
[#959]: https://github.com/YosefLab/scvi-tools/pull/959
[#966]: https://github.com/YosefLab/scvi-tools/pull/966
[#967]: https://github.com/YosefLab/scvi-tools/pull/967
[#971]: https://github.com/YosefLab/scvi-tools/pull/971
[#988]: https://github.com/YosefLab/scvi-tools/pull/988
[#989]: https://github.com/YosefLab/scvi-tools/pull/989
[#990]: https://github.com/YosefLab/scvi-tools/pull/990
[#999]: https://github.com/YosefLab/scvi-tools/pull/999
[@achille]: https://github.com/ANazaret
[@adam]: https://github.com/adamgayoso
[@adamgayoso]: https://github.com/adamgayoso
[@cane11]: https://github.com/cane11
[@casey-greene]: https://github.com/cgreene
[@cataclysmus]: https://github.com/cataclysmus
[@chenling]: https://github.com/chenlingantelope
[@david-kelley]: https://github.com/davek44
[@eddie]: https://github.com/Edouard360
[@eduardo-beltrame]: https://github.com/Munfred
[@florianbarkmann]: https://github.com/FlorianBarkmann
[@francesco-brundu]: https://github.com/fbrundu
[@gabriel]: https://github.com/gabmis
[@galen]: https://github.com/galenxing
[@galenxing]: https://github.com/galenxing
[@giovp]: https://github.com/giovp
[@gokcen-eraslan]: https://github.com/gokceneraslan
[@grst]: https://github.com/grst
[@han-yuan]: https://github.com/hy395
[@james-webber]: https://github.com/jamestwebber
[@jamie-morton]: https://github.com/mortonjt
[@jeff]: https://github.com/jeff-regier
[@jjhong922]: https://github.com/jjhong922
[@john-reid]: https://github.com/JohnReid
[@jules]: https://github.com/jules-samaran
[@marianogabitto]: https://github.com/marianogabitto
[@martinkim0]: https://github.com/martinkim0
[@max]: https://github.com/maxime1310
[@michael-raevsky]: https://github.com/raevskymichail
[@mjayasur]: https://github.com/mjayasur
[@mkarikom]: https://github.com/mkarikom
[@morris-frank]: https://github.com/morris-frank
[@munfred]: https://github.com/Munfred
[@njbernstein]: https://github.com/njbernstein
[@oscar]: https://github.com/oscarclivio
[@pierre]: https://github.com/PierreBoyeau
[@pierreboyeau]: https://github.com/PierreBoyeau
[@primoz-godec]: https://github.com/PrimozGodec
[@ricomnl]: https://github.com/ricomnl
[@rk900]: https://github.com/RK900
[@romain]: https://github.com/romain-lopez
[@romain-lopez]: https://github.com/romain-lopez
[@saketkc]: https://github.com/saketkc
[@stephen-flemming]: https://github.com/sjfleming
[@talashuach]: https://github.com/talashuach
[@tommycelsius]: https://github.com/tommycelsius
[@valentine-svensson]: https://github.com/vals
[@vitkl]: https://github.com/vitkl
[@watiss]: https://github.com/watiss
[@william-yang]: https://github.com/triyangle
[@wukathy]: https://github.com/wukathy
[@yining]: https://github.com/imyiningliu
[genomepy]: https://github.com/vanheeringen-lab/genomepy
[hatch]: https://hatch.pypa.io/latest/
[hugging face models]: https://huggingface.co/models
[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[original scbasset model]: https://github.com/calico/scBasset
[poetry]: https://python-poetry.org/
[semantic versioning]: https://semver.org/spec/v2.0.0.html
