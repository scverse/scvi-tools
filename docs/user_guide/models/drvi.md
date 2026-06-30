# DRVI

**DRVI** (Disentangled Representation Variational Inference,
Python class {class}`~scvi.external.DRVI`)
is an unsupervised deep generative model that learns an **interpretable, disentangled**
latent representation of single-cell omics data.

The advantages of DRVI are:

- Disentanglement: DRVI disentangled cell identity, signaling pathways, stress responses, developmental trajectories, and perturbation effects into distinct factors.
- Interpretability: Each factor (latent dimension) can be linked back to genes via the built-in interpretability analyses.
- Built on scVI: {class}`~scvi.external.DRVI` subclasses {class}`~scvi.model.SCVI` and only swaps
  the decoder, so it inherits the full SCVI interface — covariates, batch embedding
  (`batch_representation="embedding"`), size factors, minified mode, scArches transfer learning,
  out-of-core training via a datamodule, and the standard likelihoods (plus DRVI's log-space `pnb`).


```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/scrna/drvi`
```

```{topic} References:

-  Paper: Disentangling cellular heterogeneity into interpretable biological factors through structured latent representations.
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

Two aggregations over the split dimension are available:

- `logsumexp` (default): aggregation in log space, which makes the decoder additive in the rate.
- `mean`: average of the per-split parameters.

Two latent-to-split mappings are available (`n_latent` should be divisible by `n_split_latent`):

- `split_map` (default): the latent is reshaped into `n_split_latent` chunks and each chunk is passed through a learned per-split linear projection.
- `split_mask`: latent is split into chunks and padded with zeros so, each split `i` keeps only its own chunk.


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

## Scope: core model only

Only the **core DRVI model** and its **core interpretability methods** are included in scvi-tools
(the model, the additive decoder and log-sum-exp pooling, the log-space likelihoods, and the interpretability methods.).

The **utility functions for latent-space visualization and analysis are not** part of scvi-tools.
For those, install the full [`drvi-py`](https://github.com/theislab/drvi) package:

```bash
pip install drvi-py
```
