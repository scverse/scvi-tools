# DiagVI

**DiagVI** (Diagonal multi-modal integration Variational Inference; Python class {class}`~scvi.external.DIAGVI`) is a deep generative model for diagonal integration of unpaired multi-modal single-cell data using prior biological knowledge encoded as a guidance graph.

DiagVI is inspired by the GLUE[^ref1] architecture, which uses modality-specific variational autoencoders (VAEs) to project heterogeneous data types into a shared latent space. In contrast to GLUE’s adversarial alignment strategy, DiagVI aligns modalities using Unbalanced Optimal Transport (UOT) via the Sinkhorn divergence[^ref2], explicitly accounting for differences in cell-type composition across modalities.

The advantages of DiagVI are:

-   Flexible two-modality integration of various data types (e.g., scRNA-seq, spatial transcriptomics, spatial proteomics).
-   Full feature utilization: all features (not only overlapping ones) contribute to model training via the guidance graph and modality-specific VAEs.
-   Biologically informed alignment using prior feature correspondences via the guidance graph.
-   Robust integration of modality-specific or rare cell populations via UOT.
-   Scalable to very large datasets (>1 million cells).

The limitations of DiagVI include:

-   Currently supports integration of two modalities only.
-   Requires prior information on cross-modal feature correspondences (explicitly or implicitly).
-   May require tuning of loss weights for optimal performance (for more information, see [Practical guidance](#practical-guidance)).
-   Effectively requires a GPU for fast inference.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/multimodal/DiagVI_spatial_transcriptomics`
-   {doc}`/tutorials/notebooks/multimodal/DiagVI_spatial_proteomics`
```


## Preliminaries

### Input

DiagVI takes as input expression matrices $\mathbf{X}_1 \in \mathbb{R}^{N \times G}$ with $N$ cells and $G$ features and $\mathbf{X}_2 \in \mathbb{R}^{M \times P}$ with $M$ cells and $P$ features from two unpaired modalities.

For **count data** such as scRNA-seq data DiagVI expects a raw count expression matrix, where rows correspond to individual cells and columns correspond to features (e.g., genes).

For **continuous data** such as antibody-based single-cell proteomics data DiagVI expects a transformed (and optionally scaled) protein expression matrix, where rows correspond to individual cells and each column is a feature (e.g., a marker protein).
Input data should be preprocessed using transformations such as arcsinh or log1p, optionally followed by feature-wise scaling (e.g., z-score, min–max, or rank scaling). We recommend arcsinh transformation followed by feature-wise min–max scaling.

Optional: For both count and continuous data, DiagVI can additionally incorporate:
-   Experimental covariates such as batch annotation or confounding variables such as donor sex for one or both modalites. For simplicity, we will describe the case of categorical batch identifiers $s_n, s_m \in \{1,...,S\}$.
-   Cell label annotations that weakly inform the prior of the latent space and guide a classifier in semi-supervised training for one or both modalities. We assume categorical label identifiers $c_n, c_m \in\{1,...,C\}$.

Currently supported modalities include:
-   scRNA-seq
-   Spatial transcriptomics
-   Spatial proteomics (e.g., CITE-seq, any other antibody-based single cell measurement of protein expression)
-   Other count or continuous measurements

### Model components

DiagVI consists of several components which together define the overall training objective (see [Training Objective](#training-objective)).

#### Modality-specific variational autoencoders

DiagVI integrates two unpaired modalities by projecting their expression matrices $\mathbf{X}_1 \in \mathbb{R}^{N \times G}$ and $\mathbf{X}_2 \in \mathbb{R}^{M \times P}$ into a shared latent space in which each cell $n$ in modality 1 has a latent representation $\mathbf{z}_n \in \mathbb{R}^d$ and each cell $m$ in modality 2 has a latent representation $\mathbf{z}_m \in \mathbb{R}^d$. To find this shared state, the observed data from each modality is modeled using an independently parametrized VAE.

#### Guidance graph

Projecting cells into a common low-dimensional space alone does not result in semantically consistent embeddings in which identical cell types are assigned to the same region across both modalities. To ensure biological consistency in the latent space, a guidance graph establishes logical associations between the two modalities by connecting linked features.

Feature correspondences are encoded in the guidance graph $\mathcal{G} = (\mathcal{V}, \mathcal{E})$, where $\mathcal{V} = \mathcal{V}_1 \cup \mathcal{V}_2$ with modality-specific feature sets $\mathcal{V}_1$ and $\mathcal{V}_2$, and $\mathcal{E} \subseteq \mathcal{V} \times \mathcal{V}$ denotes the set of edges. Each edge $(i,j) \in \mathcal{E}$ is associated with

-   a weight $w_{ij} \in (0, 1]$ reflecting the confidence of the link
-   a sign $\sigma_{ij} \in \{-1,1\}$ specifying whether the interaction is associative ($\sigma_{ij} = 1$) or repressive ($\sigma_{ij} = -1$)

The graph loss encourages the inner product between embeddings $\mathbf{v}_i$ and $\mathbf{v}_j$ of linked features to be large and positive for $\sigma_{ij} = 1$ and large and negative for $\sigma_{ij} = -1$, with strength modulated by edge weights.

#### Unbalanced optimal transport

To ensure robust alignment of cells from different modalities within the shared latent space, DiagVI leverages unbalanced optimal transport (UOT). Specifically, it minimizes the de-biased Sinkhorn divergence between latent distributions using the [GeomLoss](https://www.kernel-operations.io/geomloss/#) library[^ref3].

#### Classifier

Inspired by {class}`~scvi.external.RESOLVI`[^ref4], a simple cell type classifier predicting labels $c_n$ and $c_m$ from cell latent vectors $\mathbf{z}_n$ and $\mathbf{z}_m$, respectively, is integrated into the model in the semi-supervised setting. This classifier is trained jointly with the generative model.


## Descriptive model

DiagVI assumes that observations from each modality are generated from a shared $d$-dimensional latent space. For cell $n$ and feature $g$, the the observed data $x_{ng}$ is generated conditionally on the cell latent variable $\mathbf{z}_n \in \mathbb{R}^{d}$, the feature latent variable $\mathbf{v}_g \in \mathbb{R}^{d}$, and its associated batch $s_n$.

The prior on the feature latent variable is a standard multivariate normal:

$$
\mathbf{v}_g \sim \mathcal{N}(\mathbf{0}, \mathbf{I}_d).
$$

The prior on the cell latent variable can be either a standard multivariate normal (default) or a Gaussian mixture with $L$ components:

$$
\mathbf{z}_n \sim
\begin{cases}
\mathcal{N}(\mathbf{0}, \mathbf{I}_d), \\
\sum_{l=1}^L \pi_l \, \mathcal{N}(\boldsymbol{\mu}_l, \boldsymbol{\Sigma}_l).
\end{cases}
$$

If cell type labels are provided and a Gaussian mixture prior is used,
$L$ is set to the number of unique cell types in the labeled data.


## Generative process

### Generative model

The form of the generative model depends on whether the observed data consist of discrete counts or continuous measurements. In both cases, the specific form of the data likelihood varies with the distribution used to model the observed data in each modality.

Regardless of the data type and likelihood choice, the inner product between cell and feature
latent embeddings parametrizes the decoder:

$$
\eta_{ng} = \alpha_{s_n,g} \, \mathbf{z}_n^\top \mathbf{v}_g + \beta_{s_n,g}.
$$

For count-based modalities, denoised and normalized expression proportions $\rho_{ng}$ are obtained via a softmax over features:

$$
\rho_{ng} = \frac{\exp(\eta_{ng})} {\sum_{g'} \exp(\eta_{ng'})}.
$$

The mean of the generative distribution is given by

$$
\mu_{ng} = l_n\rho_{ng},
$$

where $l_n$ denotes the observed library size of cell $n$. Counts are then modeled as

$$
x_{ng}
\sim
\mathrm{NB}(\mu_{ng}, \theta_{s_n,g})
$$

when using the negative binomial likelihood with mean $\mu_{ng}$ and dispersion $\theta_{s_n,g}$. For other likelihood choices, see [Likelihood Models](#likelihood-models).

For continuous modalities, no simplex constraint is enforced and no library size normalization is applied. When using the normal likelihood, continuous expression values are modeled as

$$
x_{ng}
\sim
\mathcal{N}(\eta_{ng}, \sigma_{s_n,g}^2).
$$

In both settings, when a batch covariate is provided, DiagVI learns batch-specific versions of scaling and bias parameters to account for batch effects within each modality. For each modality $k$ and batch $s$, the scaling and bias parameters satisfy

$$
\boldsymbol{\alpha}^{(k)}_s \in \mathbb{R}_+^{|\mathcal{V}_k|},
\qquad
\boldsymbol{\beta}^{(k)}_s \in \mathbb{R}^{|\mathcal{V}_k|}.
$$

This generative process is also summarized in the following graphical model:

:::{figure} figures/diagvi_graphical_model.svg
:align: center
:alt: DiagVI graphical model
:class: img-fluid
:width: 90%

DiagVI graphical model. Shaded nodes represent observed data, unshaded nodes represent latent variables.
:::

### Likelihood models

Depending on the modality and data characteristics, different likelihood functions can be used to reconstruct the input data. The choice of likelihood should reflect whether the data are counts or continuous measurements, and whether they exhibit excess zeros or background signal.

We generally recommend:

-   **Negative binomial (``nb``)** for count-based single-cell RNA-seq data.
-   **Normal (``normal``)** for transformed and scaled continuous measurements.

More specialized likelihoods such as zero-inflated or mixture models may be beneficial for sparse or background-contaminated measurements.

DiagVI supports the following likelihood functions:

```{eval-rst}
.. list-table::
   :widths: 10 20 30 40
   :header-rows: 1

   * - Likelihood
     - Distribution
     - Typical modality
     - Recommended preprocessing
   * - ``nb``
     - Negative Binomial
     - Count data (e.g., scRNA-seq)
     - Raw counts
   * - ``zinb``
     - Zero-Inflated Negative Binomial
     - Strongly zero-inflated count data
     - Raw counts
   * - ``nbmixture``
     - Negative Binomial Mixture
     - Protein counts with background signal (e.g., CITE-seq ADTs)
     - Raw counts
   * - ``normal``
     - Normal
     - Continuous data
     - Transformed (e.g., arcsinh or log1p) and feature-wise scaled
   * - ``log1pnormal``
     - Log1p-Normal
     - Non-negative continuous data
     - Raw data, optionally scaled
   * - ``ziln``
     - Zero-Inflated Log-Normal
     - Sparse positive continuous data
     - Raw data, optionally scaled
   * - ``zig``
     - Zero-Inflated Gamma
     - Sparse positive continuous data
     - Raw data, optionally scaled
```

### Latent variables

The exact set of variables instantiated during training depends on the chosen likelihood and on whether a Gaussian mixture prior is used for the cell latent space. Below we summarize the principal latent, decoder, and auxiliary variables used in DiagVI.

```{eval-rst}
.. list-table::
   :widths: 25 55 20
   :header-rows: 1

   * - Variable
     - Description
     - Code variable
   * - :math:`\mathbf{z}_n \in \mathbb{R}^{d}`
     - Low-dimensional latent representation of cell :math:`n`, capturing its underlying biological state.
     - ``z``
   * - :math:`\mathbf{v}_g \in \mathbb{R}^{d}`
     - Low-dimensional embedding of feature :math:`g`, inferred via the graph encoder from the guidance graph.
     - ``v``
   * - :math:`c_n \in \{1,\dots,C\}`
     - Observed cell type label of cell :math:`n` (if available), used for semi-supervised training and to inform the Gaussian mixture prior.
     - ``y``
   * - :math:`\eta_{ng} \in \mathbb{R}`
     - Decoder linear predictor obtained from the bilinear interaction between cell and feature embeddings.
     - N/A
   * - :math:`\boldsymbol{\rho}_n \in \Delta^{|\mathcal{V}_k|-1}`
     - Denoised, normalized expression proportions for cell :math:`n` (count data); constrained to the probability simplex via softmax.
     - ``px_scale``
   * - :math:`\mu_{ng} \in \mathbb{R}_+`
     - Mean of the generative distribution for feature :math:`g` in cell :math:`n`.
     - ``px_rate``
   * - :math:`\theta_{s_n,g} \in (0,\infty)`
     - Batch-specific inverse dispersion parameters used in count-based likelihoods.
     - ``px_r``
   * - :math:`\sigma_{s_n,g}^2 \in (0,\infty)`
     - Batch-specific variance parameters used in continuous likelihoods.
     - ``px_r``
   * - :math:`\delta_g \in (0,\infty)`
     - Dropout or zero-inflation parameters (used for zero-inflated likelihoods).
     - ``px_dropout``
   * - :math:`\pi_l`
     - Mixture weights of the Gaussian mixture prior (if enabled).
     - ``gmm_logits``
```


## Inference

Since the posterior distributions over the latent variables are intractable, DiagVI uses variational inference and specifically auto-encoding variational Bayes (see {doc}`/user_guide/background/variational_inference`) to jointly learn model parameters and approximate posterior distributions.

For a given modality, the variational distribution factorizes as

$$
q_\eta(\mathbf{z}_n, \mathbf{v} \mid x_n, \mathcal{G}) = q_\eta(\mathbf{z}_n \mid x_n) \; q_\eta(\mathbf{v} \mid \mathcal{G}).
$$

Here $\eta$ is a set of parameters corresponding to inference neural networks (encoders).


## Training objective

DiagVI is trained by minimizing a weighted sum of loss terms (weighted by `lam_*` parameters) corresponding to the [Model Components](#model-components) introduced above:

-   **Modality-specific VAEs**: The data reconstruction loss (`lam_data`) measures how well each modality-specific decoder reconstructs the observed data, while the KL divergence (`lam_kl`) regularizes the cell latent variables by encouraging adherence to the prior.
-   **Guidance graph**: The graph reconstruction loss (`lam_graph`) enforces biological consistency between feature embeddings using the guidance graph.
-   **UOT**: The UOT alignment loss (`lam_sinkhorn`) aligns cell distributions across modalities using unbalanced optimal transport.
-   **Classifier**: The classification loss (`lam_class`, optional) enables supervised or semi-supervised training via cell-type labels.

The `lam_*` parameters control the relative importance of within-modality reconstruction (`lam_data`, `lam_kl`), cross-modality alignment (`lam_graph`, `lam_sinkhorn`), and optional supervision (`lam_class`). DiagVI uses sensible defaults for all of these values, but they may require further tuning depending on the type of data being integrated (see [Practical guidance](#practical-guidance)).


## Practical guidance

### Loss weights

DiagVI provides default values for all loss weights that have been tested across multiple datasets and integration tasks. However, depending on the modalities and the integration setting, tuning some of these weights can improve performance.

In practice, we recommend paying particular attention to:

-   `lam_sinkhorn` (cross-modality alignment strength), and
-   `lam_class` (classification strength, if cell-type labels are provided).

These two weights control the trade-off between alignment across modalities and separation of cell types within the latent space.

When integrating **very different modalities** (e.g., scRNA-seq + spatial proteomics), stronger alignment is typically required. In such cases:

-   Higher `lam_sinkhorn`
-   Lower `lam_class` (if labels are used)

This encourages the model to prioritize matching global cell distributions across modalities, even if the feature spaces differ substantially.

When integrating **more similar modalities** (e.g., scRNA-seq + spatial transcriptomics), less alignment pressure is needed. In such cases:

-   Lower `lam_sinkhorn`
-   Higher `lam_class` (if labels are used)

Here, the modalities are already structurally similar, so stronger supervision can help maintain clearer separation of biologically distinct cell types.

As a general rule, we recommend `lam_sinkhorn` and `lam_class` values between 1 and 100, which work well in most settings.

The default loss weight values are:

-   `lam_data` = 0.1
-   `lam_kl` = 1.0
-   `lam_graph` = 15.0
-   `lam_sinkhorn` = 85.0
-   `lam_class` = 5.0

These defaults place relatively strong emphasis on cross-modality alignment and are therefore well suited for integrating heterogeneous modalities. For more similar modalities, we recommend reducing `lam_sinkhorn` (e.g., 5-20) and increasing `lam_class` (e.g., 60-80), depending on label availability and desired separation.

In practice, small grid searches over `lam_sinkhorn` and `lam_class` are often sufficient to find well-performing settings.

### Sinkhorn parameters

Three parameters determine the behavior of the Sinkhorn divergence used for cross-modality alignment:

-   `p`: order of the ground cost (Wasserstein-p distance). The default is p = 2, corresponding to a squared Euclidean ground cost.
-   `blur`: controls the strength of entropic regularization and therefore the smoothness of the transport plan.
-   `reach`: controls the marginal relaxation parameter in unbalanced optimal transport, penalizing deviations from strict mass conservation.

Although these parameters can be specified manually, we generally **do not recommend modifying them**.

By default, DiagVI follows the heuristic strategy introduced in [OTT-JAX](https://ott-jax.readthedocs.io/). In this setting:
-   `blur` is computed adaptively at each optimization step from the minibatch cost matrix, and
-   `reach` is derived as a function of blur.

This adaptive strategy makes the Sinkhorn loss scale-aware and robust across datasets, and in practice removes the need for manual tuning. Adjusting these parameters is rarely necessary.


## Tasks

Here we provide an overview of some of the tasks that DiagVI can perform. Please see {class}`~scvi.external.DIAGVI` for the full API reference.

### Dimensionality reduction

```
>>> latents = model.get_latent_representation()
>>> adata_rna.obsm["X_diagvi"] = latents["rna"]
>>> adata_protein.obsm["X_diagvi"] = latents["protein"]
```
Aligned representations can be used jointly for downstream analysis:

```
>>> import scanpy as sc
>>> # Use RNA latent for clustering
>>> sc.pp.neighbors(adata_rna, use_rep="X_diagvi")
>>> sc.tl.umap(adata_rna)
>>> sc.tl.leiden(adata_rna)
```

### Cross-modal feature imputation

DiagVI can impute features from one modality to another using {func}`~scvi.external.DIAGVI.get_imputed_values`.

```
>>> # Impute protein expression from RNA data
>>> imputed_protein = model.get_imputed_values(
...     source_name="rna",
...     source_adata=adata_rna,
... )
>>> adata_rna.obsm["imputed_protein"] = imputed_protein
```

You can also specify target batch and library size for counterfactual predictions:

```
>>> # Impute with specific target batch
>>> imputed = model.get_imputed_values(
...     source_name="rna",
...     target_batch="batch_1",
...     target_libsize=10000,
... )
```

### Cell label transfer


## References

[^ref1]:
    Cao, Zhi-Jie; Gao, Ge (2022),
    _Multi-omics single-cell data integration and regulatory inference with graph-linked embedding_,
    [Nature Biotechnology](https://www.nature.com/articles/s41587-022-01284-4).

[^ref2]:
    Séjourné, Thibault; Feydy, Jean; Vialard, François-Xavier; Trouvé, Alain; Peyré, Gabriel (2023),
    _Sinkhorn Divergences for Unbalanced Optimal Transport_,
    [arXiv](https://arxiv.org/abs/1910.12958).

[^ref3]:
    Feydy, Jean; Séjourné, Thibault; Vialard, François‑Xavier; Amari, Shun‑ichi; Trouvé, Alain; Peyré, Gabriel (2019)
    _Interpolating between Optimal Transport and MMD using Sinkhorn Divergences_
    [The 22nd International Conference on Artificial Intelligence and Statistics](https://arxiv.org/abs/1810.08278).

[^ref4]:
    Adam Gayoso\*, Zoë Steier\*, Romain Lopez, Jeffrey Regier, Kristopher L Nazor, Aaron Streets, Nir Yosef (2021),
    _Joint probabilistic modeling of single-cell multi-omic data with totalVI_,
    [Nature Methods](https://www.nature.com/articles/s41592-020-01050-x).
