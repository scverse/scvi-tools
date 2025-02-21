# Frequently asked questions

This section contains answers to frequently asked questions. Don't see what you're looking for
here? Feel free to check out our [Discourse](https://discourse.scverse.org/c/help/scvi-tools/)
or chat with us directly on [Zulip](https://scverse.zulipchat.com/#narrow/stream/324229-scvi-tools).

## What is the difference between `batch_key` and `categorical_covariate_keys`?

While both `batch_key` and `categorical_covariate_keys` are meant to be used for categorical
covariates whose effects we want to regress out, they have support for different features.
`batch_key` generally supports more features and should be used before `categorical_covariate_keys`
if possible. In particular, `batch_key` supports the following compared to
`categorical_covariate_keys`:

- Used to specify the main technical effect in the data, such as the sequencing lab, lane,
    sequencing type, or dataset of origin, for example.
- Support for per-category dispersion parameters by specifying `dispersion="gene-batch"` for
    models such as [`SCVI`](https://docs.scvi-tools.org/en/latest/api/reference/scvi.model.SCVI.html#scvi.model.SCVI).
- Support for flexible embedding modeling by specifying `batch_representation="embedding"`.
- Support for counterfactual decoding by specifying the `transform_batch` argument in methods such
    as [`get_normalized_expression`](https://docs.scvi-tools.org/en/latest/api/reference/scvi.model.SCVI.html#scvi.model.SCVI.get_normalized_expression).
- For models where `use_observed_lib_size=False` can be set, the model will treat the library size
    as a learned parameter, where the variational prior is set to a normal distribution with mean
    and variance computed per-batch.

On the other hand, `categorical_covariate_keys` supports the following features compared to
`batch_key`:

- Definition of multiple covariates at once (e.g.,
`categorical_covariate_keys=["assay_type", "donor"]`).

They also share some features:

- One-hot encoded by default
- Only passed into the decoder by default, and have the same behavior when specifying
    `encode_covariates` and `deeply_inject_covariates`.
- Meant to capture technical nuisance effects.

## My model errors out during training due to `NaN`s - how can I fix this?

There is no straightforward solution since these are usually caused by numerical instabilities, but
here are several factors to consider:

- **Data quality**: Ensure that the data is appropriately preprocessed (_e.g._ removing cells/spots
    with low counts).
- **Data distribution**: Ensure that, if the model expects raw counts (_e.g._ scVI), the data is
    not normalized or transformed in any way. This often means passing in the correct `layer`
    argument to the `setup_anndata` method.
- **Training parameters**: You may need to adjust parameters such as the learning rate or batch
    size to stabilize training, or use gradient clipping to prevent exploding gradients.
- **Model activations**: Some models contain exponential activations due to reproducibility
    reasons, which may be sometimes numrically unstable. Consider using another activation such as
    the softplus for more stability.
- **Adversarial training**: For models with adversarial training (_e.g._ totalVI), consider turning
    off the component to see if the issue is resolved.
- **Using SaveCheckpoint Callback**: Starting v1.3.0, we added the on_exception option to the callback, {class}`scvi.train._callbacks.SaveCheckpoint`, that in case of model error exception during training, resulted from Nan's in loss or gradients, will save the best model ("best" in terms of what is the monitored metric) that was achieved upto the point of failure. It will not gracefully shutdown, but the user responsibility will be to load it back and continue the analysis, e.g., user can take the optimal model that was saved, or continue training it with perhaps different parameters, where, in such case it might continue to train passed the point of previous failure. Note this option will add some overhead to train time. It can be used by adding the following parameter to the train function: `callbacks=[SaveCheckpoint(monitor="elbo_validation", load_best_on_end=True)]`
