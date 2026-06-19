# DRVI

**DRVI** (Disentangled Representation Variational Inference,
Python class {class}`~scvi.external.DRVI`)
is an unsupervised deep generative model that learns an **interpretable, disentangled**
latent representation of single-cell omics data.

The advantages of DRVI are:

- Interpretable factors: Each latent dimension tends to capture an independent, biologically
  meaningful axis of variation, which can be linked back to genes via the built-in
  interpretability analyses.
- Built on scVI: {class}`~scvi.external.DRVI` subclasses {class}`~scvi.model.SCVI` and only swaps
  the decoder, so it inherits the full SCVI interface — covariates, batch embedding
  (`batch_representation="embedding"`), size factors, minified mode, scArches transfer learning,
  out-of-core training via a datamodule, and the standard likelihoods (plus DRVI's log-space `pnb`).

The limitations of DRVI include:

- Interpretability assumes one split per latent dimension (`n_split_latent == n_latent`, the
  default).

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/scrna/drvi`
```

```{topic} References:

-  Paper: Moinfar and Theis. Unsupervised deep disentangled representation of single-cell omics.
bioRxiv (2024): https://doi.org/10.1101/2024.11.06.622266
```

## Method background

DRVI is a variational autoencoder (VAE): the integrated, disentangled representation corresponds to
the latent-space embedding of the cells. It shares the encoder, prior (`N(0, I)`) and reconstruction
likelihoods with {class}`~scvi.module.VAE`. Disentanglement is induced **in the decoder** rather
than through the prior.

### Disentanglement via an additive split decoder

The latent vector of dimension `K` is mapped into `n_split_latent` independent groups ("splits").
Each split is decoded **independently** and the per-split outputs are aggregated to form the final
reconstruction parameters. Because every split can only influence the output through this shared
aggregation, the model is encouraged to allocate independent factors of variation to different
latent dimensions.

Two latent-to-split mappings are available:

- `split_diag`: requires `n_split_latent == n_latent`; split `i` sees only latent dimension `i`
  (via a diagonal embedding).
- `split_map` (default): the latent is reshaped into `n_split_latent` chunks and each chunk is
  passed through a learned per-split linear projection
  ({class}`stacked_linear.StackedLinearLayer`).

Two aggregations over the split dimension are available:

- `logsumexp` (default): aggregation in log space, which makes the decoder additive in the rate.
- `mean`: average of the per-split parameters.

### Likelihoods modeled in log space

DRVI models the generative parameters in **log space**, which is what makes the `logsumexp`
aggregation a genuine additive decoder (summing the per-split contributions of the rate). In
addition to scvi's `nb`, `zinb` and `poisson`, DRVI adds:

- `pnb` — *parametrized negative binomial*: the same mean as `nb` (``library * softmax``) but
  parametrized and evaluated in log space via {class}`~scvi.external.drvi.LogNegativeBinomial`,
  which is numerically stable and composes the per-split log contributions directly.
- `normal` / `normal_unit_var` — a Gaussian whose mean is modeled **directly** (no library/softmax),
  with the per-gene variance modeled in log space (fixed to one for `normal_unit_var`); intended for
  continuous, e.g. log-normalized, input.

### Interpretability

Because the decoder is additive over splits, the contribution of each split to the reconstructed
expression of every gene can be measured directly. DRVI exposes per-split parameters in an
inspection mode and provides methods to quantify these contributions both in-distribution (over the
observed cells) and out-of-distribution (by traversing each latent dimension), and to extract the
top genes per latent dimension. See {meth}`~scvi.external.DRVI.set_latent_dimension_stats`,
{meth}`~scvi.external.DRVI.calculate_interpretability_scores` and
{meth}`~scvi.external.DRVI.get_interpretability_scores`.
