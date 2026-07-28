# gimVI

**gimVI** {cite:p}`Lopez19` (Generative Imputation Model with Variational Inference; Python class {class}`~scvi.external.GIMVI`) is a deep generative model for the joint analysis of an scRNA-seq dataset and a spatial transcriptomics dataset that measures only a subset of the genes. gimVI learns a single shared latent space for the two modalities and uses it to **impute the genes that are missing from the spatial assay** by transferring information from the (genome-wide) sequencing reference.

The two datasets are modeled by a single joint variational autoencoder ({class}`~scvi.external.gimvi.JVAE`) with modality-specific encoder and decoder heads but a shared latent representation. Because the spatial panel is a subset of the genes captured by scRNA-seq, the decoder reconstructs into the full sequencing gene space, and each modality's reconstruction loss is masked to only the genes it actually observes. Alignment of the two modalities in the shared latent space is encouraged by an adversarial classifier that is trained to predict the modality of origin from a cell's latent code, while the generative model is trained to fool it.

The advantages of gimVI are:

-   Imputation of genes that are not measured in the spatial assay, using a paired (genome-wide) scRNA-seq reference.
-   A shared latent space that integrates spatial and sequencing cells for joint visualization and clustering.
-   Flexible, per-modality generative distributions (e.g. `zinb` for the sequencing data and `nb` for the spatial data) and per-modality choice of whether to model library size.
-   Handles the two assays measuring different gene sets and different numbers of cells.

The limitations of gimVI include:

-   Supports exactly two datasets (one sequencing, one spatial).
-   The spatial genes must be a subset of the sequencing genes.
-   Imputation quality depends on how well the scRNA-seq reference matches the spatial tissue.

:::{note}
Starting scVI-Tools v1.5 this model is part of scVIVA-Tools, and no longer being maintained here.
:::

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/spatial/gimvi_tutorial`
```

:::{note}
gimVI is a spatial transcriptomics model that will be moved to the scvi-tools spatial companion package `scviva-tools` starting in scvi-tools v1.5 and will no longer be supported in scvi-tools; it is scheduled for deprecation in v1.6.
:::


## Preliminaries

gimVI takes as input two raw-count gene expression matrices:

-   a sequencing matrix $X^{(0)} \in \mathbb{N}^{N_0 \times G_0}$ with $N_0$ cells and $G_0$ genes (the full gene set), and
-   a spatial matrix $X^{(1)} \in \mathbb{N}^{N_1 \times G_1}$ with $N_1$ cells and $G_1$ genes,

where the spatial genes are a subset of the sequencing genes ($\mathcal{G}_1 \subseteq \mathcal{G}_0$). Internally the decoder operates on the full sequencing gene space of size $G = G_0$, and a per-modality index mapping $\mathcal{M}_m$ records where each modality's observed genes live in that space. The spatial reconstruction is therefore evaluated only on the columns corresponding to $\mathcal{G}_1$.

Each dataset can additionally carry a categorical batch covariate $s$ and, optionally, categorical cell-type labels $c$. Both datasets are registered independently with two separate calls to {meth}`~scvi.external.GIMVI.setup_anndata` before the joint model is constructed:

```
>>> scvi.external.GIMVI.setup_anndata(adata_seq, batch_key="batch", labels_key="labels")
>>> scvi.external.GIMVI.setup_anndata(adata_spatial, batch_key="batch", labels_key="labels")
>>> model = scvi.external.GIMVI(adata_seq, adata_spatial)
>>> model.train()
```

The two datasets share a single batch space: batch indices of the spatial dataset are offset by the number of sequencing batches so that the two modalities never collide on a batch id (see [#2446](https://github.com/scverse/scvi-tools/issues/2446)).

### Model architecture

gimVI is built from the following components, which together define the [training objective](#training-objective):

#### Joint encoder

A single multi-head encoder ({class}`~scvi.nn.MultiEncoder`) maps both modalities into a shared $d$-dimensional latent space. Each modality has its own input layer(s) (because the two assays have different numbers of genes), followed by shared layers that produce the variational posterior $q(\mathbf{z}_n \mid x_n)$. The `mode` argument selects which input head is used.

For each modality, an optional, separate library-size encoder produces a latent library size. Whether a modality models library size as a latent variable or instead uses the observed log total counts is controlled by `model_library_size` (default `[True, False]`: the sequencing data models library size, the spatial data uses the observed value).

#### Joint decoder

A single, shared multi-decoder ({class}`~scvi.nn.MultiDecoder`) maps the shared latent code back into the full sequencing gene space of size $G$, conditioned on the batch covariate (the decoder weights are **not** modality-specific). Because the decoder reconstructs all $G$ genes for both modalities, the latent representation of a spatial cell can be decoded into genes that were never measured spatially — this is exactly the imputation mechanism (see [Imputation of missing spatial genes](#imputation-of-missing-spatial-genes)). The modality index only selects which gene-index mapping is used to renormalize the decoded expression over each modality's observed genes; it does not route through different decoder networks.

#### Adversarial modality classifier

A small classifier ({class}`~scvi.module.Classifier`, 3 layers, 2 output classes) is trained to predict, from a latent code $\mathbf{z}$, which modality it came from. The generative model is trained adversarially to *fool* this classifier, which pushes the sequencing and spatial latent distributions to overlap and yields an integrated, modality-mixed latent space. This is implemented in {class}`~scvi.external.gimvi._task.GIMVITrainingPlan` (a subclass of {class}`~scvi.train.AdversarialTrainingPlan`), and the strength of the adversarial term is set by the `kappa` argument of {meth}`~scvi.external.GIMVI.train`.


## Generative process

For a cell $n$ in modality $m \in \{0, 1\}$ (0 = sequencing, 1 = spatial) with batch $s_n$, gimVI posits the following generative process over the full gene space:

```{math}
:nowrap: true

\begin{align}
 \mathbf{z}_n &\sim \mathrm{Normal}(0, I) \\
 \ell_n &\sim
 \begin{cases}
   \mathrm{LogNormal}(\ell_\mu^\top s_n,\ \ell_{\sigma^2}^\top s_n) & \text{if modality } m \text{ models library size} \\
   \delta(\log \textstyle\sum_g x_{ng}) & \text{otherwise (observed library size)}
 \end{cases} \\
 \rho_n &= f_w(\mathbf{z}_n, s_n) \\
 x_{ng} &\sim \mathrm{ObservationModel}_m(\ell_n \rho_{ng},\ \theta_{s_n g}),\quad g \in \mathcal{G}_m
\end{align}
```

The latent code $\mathbf{z}_n$ is **shared across modalities** — the same prior and the same latent space are used for both datasets. The shared decoder $f_w$ produces normalized expression $\rho_n$ over the full gene set; for modality $m$, $\rho_n$ is renormalized to sum to 1 over only the genes observed in that modality ($\mathcal{M}_m$), and the mean rate is $\mu_{ng} = \ell_n \rho_{ng}$.

The observation model is chosen per modality via `generative_distributions` (default `["zinb", "nb"]`). gimVI supports:

-   `zinb` — Zero-Inflated Negative Binomial (default for the sequencing data),
-   `nb` — Negative Binomial (default for the spatial data),
-   `poisson` — Poisson.

The inverse dispersion $\theta$ is by default gene- and batch-specific (`dispersion="gene-batch"`); `gene` and `gene-label` parameterizations are also available.

The latent variables are summarized below:

```{eval-rst}
.. list-table::
   :widths: 25 60 15
   :header-rows: 1

   * - Latent variable
     - Description
     - Code variable
   * - :math:`\mathbf{z}_n \in \mathbb{R}^{d}`
     - Shared low-dimensional representation of a cell, common to both modalities.
     - ``z``
   * - :math:`\rho_n \in \Delta^{|\mathcal{M}_m|-1}`
     - Denoised, normalized expression for cell :math:`n`, renormalized over the genes observed in its modality.
     - ``px_scale``
   * - :math:`\ell_n \in (0, \infty)`
     - Library size. Modeled as a latent variable for modalities flagged in ``model_library_size``; otherwise set to the observed log total count.
     - ``library``
   * - :math:`\mu_{ng} \in (0, \infty)`
     - Mean of the count likelihood, :math:`\ell_n \rho_{ng}`.
     - ``px_rate``
   * - :math:`\theta_{s_n g} \in (0, \infty)`
     - Inverse dispersion of the (zero-inflated) negative binomial, by default gene- and batch-specific.
     - ``px_r``
```


## Inference

gimVI uses variational inference (auto-encoding variational Bayes; see {doc}`/user_guide/background/variational_inference`) with the per-modality factorization

```{math}
:nowrap: true

\begin{align}
   q_\eta(\mathbf{z}_n, \ell_n \mid x_n) := q_\eta(\mathbf{z}_n \mid x_n)\, q_\eta(\ell_n \mid x_n).
\end{align}
```

The latent posterior $q_\eta(\mathbf{z}_n \mid x_n)$ is produced by the shared multi-head encoder (using the input head for the cell's modality), and $q_\eta(\ell_n \mid x_n)$ by the modality's optional library encoder. Inputs are log-transformed (`log(1 + x)`) before encoding for numerical stability.


## Training objective

gimVI minimizes a per-modality evidence lower bound — a reconstruction term plus KL terms on the latent variables — summed across both modalities, with an adversarial term added on top:

-   **Reconstruction**: the negative log-likelihood of the observed counts under the chosen observation model, **masked to the genes observed in each modality** (the spatial term ignores genes outside the spatial panel).
-   **KL divergence on $\mathbf{z}$**: regularizes the shared latent variable toward the standard normal prior.
-   **KL divergence on $\ell$**: applied only for modalities that model library size as a latent variable, with a batch-specific log-normal prior.
-   **Adversarial loss**: the modality classifier is trained to predict modality from $\mathbf{z}$, while the generative model is trained to fool it. The relative weight is controlled by `kappa` (the `scale_adversarial_loss` of the training plan); when set to `"auto"` it follows `1 - kl_weight` over the warmup schedule.

The two objectives (model vs. classifier) are optimized in alternating steps with manual optimization, as in {class}`~scvi.train.AdversarialTrainingPlan`.

### Cyclic multi-dataset data loading

The two datasets generally have different numbers of cells, so a custom data loader, {class}`~scvi.external.gimvi._task.CyclicMultiDataLoader`, combines their per-modality dataloaders: it iterates the longest loader once per epoch and **cycles the shorter loader(s)** so that every training step yields one minibatch from each modality. It accepts either a sequence (yielding a tuple of minibatches) or a mapping of modality name → loader (yielding a dict, preserving modality names).

This loader is shared infrastructure: {class}`~scvi.external.DIAGVI` imports it (`from scvi.external.gimvi._task import CyclicMultiDataLoader as TrainDL`) and uses the mapping form to drive its own multi-modality training, so the same cycling logic backs both gimVI and DiagVI.


## Tasks

Here we provide an overview of some of the tasks that gimVI can perform. Please see {class}`~scvi.external.GIMVI` for the full API reference.

### Dimensionality reduction

{meth}`~scvi.external.GIMVI.get_latent_representation` returns the shared latent representation for **each** dataset as a list `[latent_seq, latent_spatial]`. By default the mean of the approximate posterior is returned (`deterministic=True`).

```
>>> latent_seq, latent_spatial = model.get_latent_representation()
>>> adata_seq.obsm["X_gimvi"] = latent_seq
>>> adata_spatial.obsm["X_gimvi"] = latent_spatial
```

Because the latent space is shared and modality-mixed (via the adversarial classifier), the two representations can be concatenated for a joint neighbor graph and UMAP to assess integration quality:

```
>>> import scanpy as sc, anndata as ad, numpy as np
>>> joint = ad.AnnData(np.concatenate([latent_seq, latent_spatial]))
>>> joint.obs["modality"] = ["seq"] * len(latent_seq) + ["spatial"] * len(latent_spatial)
>>> sc.pp.neighbors(joint, use_rep="X")
>>> sc.tl.umap(joint)
```

### Imputation of missing spatial genes

The central task of gimVI is to impute, for the spatial cells, the genes that were only measured in the scRNA-seq data. {meth}`~scvi.external.GIMVI.get_imputed_values` returns the decoded expression for **all** sequencing genes for each dataset:

```
>>> imputed_seq, imputed_spatial = model.get_imputed_values(normalized=True)
>>> # imputed_spatial now contains values for every sequencing gene,
>>> # including those absent from the original spatial panel.
```

-   `normalized=True` (default) returns the normalized expression frequencies $\rho_n$ (decoder `sample_scale`); `normalized=False` returns the unnormalized mean rate $\mu_n = \ell_n \rho_n$ (decoder `sample_rate`).
-   `decode_mode` selects which dataset's decoding is applied. The cell is always encoded with its own modality's encoder; gimVI then uses a **single shared decoder**, and `decode_mode` only changes which modality's gene-index mapping is used to renormalize the decoded expression over the probability simplex (`indices_mappings[decode_mode]`). It does not switch to a separate per-modality decoder network.

A common evaluation (used in the tutorial) holds out a fraction of genes from the spatial panel, imputes them with gimVI, and measures the per-gene Spearman correlation between the held-out true spatial expression and the imputed values.

### Saving and loading

gimVI uses a custom {meth}`~scvi.external.GIMVI.save` / {meth}`~scvi.external.GIMVI.load` that persists the two registries and both sets of `var_names` (the model stores two AnnData managers). `save_anndata=True` writes both `adata_seq` and `adata_spatial` alongside the model. Legacy saves (pre-v0.15.0) can be upgraded with {meth}`~scvi.external.GIMVI.convert_legacy_save`.
