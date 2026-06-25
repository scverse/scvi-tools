# VeloVI

**VeloVI** {cite:p}`GayosoWeiler23` (Python class {class}`~scvi.external.VELOVI`) is a
variational model for RNA velocity. It models paired spliced and unspliced RNA abundance
to infer latent time, gene-specific kinetic rates, cell-level velocities, expression
fits, and velocity uncertainty.

The advantages of VeloVI are:

-   It provides a probabilistic RNA velocity model with posterior sampling.
-   It returns velocity estimates, latent time, state assignments, fitted expression, and
    learned kinetic rates through the scvi-tools API.
-   It includes utilities for intrinsic directional uncertainty and permutation scores.

The limitations of VeloVI include:

-   It requires preprocessed spliced and unspliced layers, usually from an RNA velocity
    workflow such as scVelo preprocessing.
-   The kinetic interpretation depends on the adequacy of the assumed induction and
    repression dynamics for the selected genes.
-   It is not a count-based scRNA-seq expression model and is not a replacement for
    general-purpose models such as scVI.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/scrna/velovi`
```

## Preliminaries

VeloVI takes as input an AnnData object with two registered layers:

-   a spliced abundance layer, passed as `spliced_layer`, and
-   an unspliced abundance layer, passed as `unspliced_layer`.

These layers are registered with {meth}`~scvi.external.VELOVI.setup_anndata` and are
treated as continuous preprocessed values rather than raw count data:

```
>>> VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
>>> model = VELOVI(adata)
>>> model.train()
```

The tutorial uses scVelo-style moment-smoothed layers (`"Ms"` and `"Mu"`), then stores
the model outputs back into `adata.layers`, `adata.var`, and `adata.obs`.

## Model Overview

VeloVI encodes the concatenated spliced and unspliced abundance vectors for each cell into
a latent representation $z_i$. A decoder maps $z_i$ to gene- and cell-specific parameters
for RNA velocity dynamics.

For each gene, VeloVI models induction and repression dynamics with transcription,
splicing, and degradation rates:

```{math}
:nowrap: true

\begin{align}
 \frac{du}{dt} &= \alpha(t) - \beta u \\
 \frac{ds}{dt} &= \beta u - \gamma s,
\end{align}
```

where $u$ is unspliced abundance, $s$ is spliced abundance, $\alpha$ is the transcription
rate, $\beta$ is the splicing rate, and $\gamma$ is the degradation rate.

The decoder outputs mixture weights over four kinetic states:

-   induction,
-   induction steady state,
-   repression, and
-   repression steady state.

It also outputs time parameters that are used with gene-specific switch times and kinetic
rates to compute the expected spliced and unspliced abundances.

## Inference

VeloVI uses amortized variational inference. The loss combines reconstruction terms for
spliced and unspliced abundance, a KL penalty for the latent representation, a KL penalty
for the state-mixture distribution, and a penalty that encourages the induction endpoint
to match the learned switch state.

The learned rates are returned by {meth}`~scvi.external.VELOVI.get_rates`, which includes
the splicing rate `beta`, degradation rate `gamma`, transcription on-state rate `alpha`,
transcription off-state rate `alpha_1`, and switching rate `lambda_alpha`.

## Tasks

Here we provide an overview of common tasks. Please see {class}`~scvi.external.VELOVI`
for the full API reference.

### Latent Representation

The VeloVI latent representation can be used for visualization and downstream
neighborhood analyses:

```
>>> adata.obsm["X_velovi"] = model.get_latent_representation()
```

### Latent Time and Velocity

{meth}`~scvi.external.VELOVI.get_latent_time` returns cell-by-gene latent time estimates,
and {meth}`~scvi.external.VELOVI.get_velocity` returns velocity estimates. Both methods
support posterior sampling and gene subsetting.

```
>>> latent_time = model.get_latent_time(n_samples=25)
>>> velocity = model.get_velocity(n_samples=25, velo_statistic="mean")
```

### State Assignment and Expression Fit

{meth}`~scvi.external.VELOVI.get_state_assignment` returns probabilities over the four
kinetic states for each cell and gene, with an option for hard state assignments.
{meth}`~scvi.external.VELOVI.get_expression_fit` returns fitted spliced and unspliced
abundances from the model.

### Uncertainty and Permutation Scores

{meth}`~scvi.external.VELOVI.get_directional_uncertainty` samples velocities from the
posterior and computes cell-level uncertainty in velocity direction. The
{meth}`~scvi.external.VELOVI.get_permutation_scores` method compares expression fits
against data with spliced and unspliced layers shuffled within cell-type labels, producing
gene-by-cell-type dynamical scores.
