# CytoVI

**CytoVI** [^ref1] (Python class {class}`~scvi.external.CYTOVI`) is a generative model for cytometry data that leverages deep probabilistic latent variable modeling to enable denoising, imputation, integration, and differential analysis across technologies and batches.

The advantages of CytoVI are:

- Provides a denoised and batch-corrected low-dimensional cell state representations across antibody-based single cell technologies.
- Facilitates the integration of data from different antibody panels and the imputation of non-overlapping markers.
- Enables cross-technology integration (e.g. flow cytometry, mass cytometry and CITE-seq) and modality imputation, facilitating unified analysis across cytometry and transcriptomics.
- Allows for the automated identification of disease-associated cell states (label-free differential abundance analysis) and differential protein expression between groups of cells.
- Scalable to very large datasets (>20 million cells).

The limitations of CytoVI include:

- Requires at least partial feature overlap across datasets for effective integration.
- Effectively requires a GPU for fast inference.
- Assumes measurements have been corrected for fluorescent spillover and preprocessed according to standard practice.

```{topic} Tutorials:
- {doc}`/tutorials/notebooks/cytometry/CytoVI_batch_correction_tutorial`
- {doc}`/tutorials/notebooks/cytometry/CytoVI_advanced_tutorial`
```

<!--  note: to populate -->

## Preliminaries
CytoVI can process protein expression matrices from the following technologies:
- Flow cytometry (traditional PMT-based cytometers or full spectrum analyzers)
- Mass cytometry
- CITE-seq (providing RNA and protein expression for the same cells)
- Any other antibody-based single cell measurement of protein expression


For each of these technologies CytoVI expects as input:
- A transformed (and optionally scaled) protein expression matrix $X \in \mathbb{R}^{N \times P}$, where each row is a single cell among $N$ total cells and each column is a protein marker among $P$ total proteins.
- Optionally, experimental covariates that the user wishes to control for such as batch annotations, different technologies or confounding variables such as donor sex. For simplicity we will describe the case of a one-hot encoded batch identifier $s_n \in \{1,...,S\}$.
- Optionally, cell label annotations can be utilized to weakly inform the prior of the latent space. We assume a one-hot encoded label identifier $y_n \in\{1,...,Y\}$.
- Optionally, sample annotations (e.g. indicating which cells were measured from which patient) to inform differential abundance and expression analyses.


Preprocessing expected:
CytoVI expects the input matrix $x_{np}$ to be processed using common preprocessing options for cytometry data in the field. These include for instance arcsinh, log1p, biexponential or logicle transformations and are optionally followed by feature-wise scaling (e.g. z-score, min-max or rank-scaled).


## Descriptive model
CytoVI is a latent variable model that assumes each cell $n$ has a $d$-dimensional latent representation $z_n$ capturing its intrinsic state that is decoupled from the variation in $s$. It models the observed protein expression $x_{np}$ for protein $p$ as a function of $z_n$ and its associated batch $s_n$.

We assume a prior on the latent space $z_n$ that can be either:

$$
z_n \sim
\begin{cases}
\mathcal{N}(0, 1), & \text{if isotropic Gaussian prior} \\
\sum_{k=1}^K \pi_k \, \mathcal{N}(\mu_k, \sigma_k^2), & \text{if mixture of Gaussians prior}
\end{cases}
$$


<!-- ```{math}
:nowrap: true
\begin{equation}
z_n \sim
\begin{cases}
\mathcal{N}(0, 1), & \text{if isotropic Gaussian prior} \\
\sum_{k=1}^K \pi_k \, \mathcal{N}(\mu_k, \sigma_k^2), & \text{if mixture of Gaussians prior}
\end{cases}
\end{equation}
``` -->

<!-- to do: exchange all math to the above syntax -->

By default, a mixture of Gaussians is used, which allows more expressive modeling and improves integration of heterogeneous data. Optionally, prior weights can be informed by known labels $y$:

$$
\pi_k' = \pi_k + \lambda y
$$

where $\lambda$ is a tunable weight (default $\lambda=10$).


## Generative process

The decoder network maps latent variables and batch vectors to a vector of protein expression parameters:

$$
f_x: \mathbb{R}^d \times \{0,1\}^S \to \mathbb{R}^P
$$

We assume the observed expression $x_{np}$ follows one of the following distributions:

$$
x_{np} \mid z_n, s_n \sim
\begin{cases}
\mathcal{N}(\mu_{np}, \sigma_{np}^2), & \text{if Gaussian model (default)} \\
\text{Beta}(\alpha_{np}, \beta_{np}), & \text{if Beta model}
\end{cases}
$$

The Gaussian likelihood is suited for arcsinh/log-transformed (and optionally scaled) cytometry data. The Beta likelihood requires a scaling of the data to a $[0, 1]$ range.

## Inference
Since the posterior $p(z_n \mid x_n)$ is intractable, we use variational inference and specifically auto-encoding variational bayes (see {doc}`/user_guide/background/variational_inference`) to learn both the model parameters (the
neural network params, params for the protein observation model, etc.) and an approximate posterior distribution.

## Handling of overlapping antibody panels
To integrate datasets with different protein panels, CytoVI employs a masking strategy inspired by {class}`~scvi.model.TOTALVI` [^ref2].

Let $x_n^{(s)} \in \mathbb{R}^{P_s}$ be the observed input vector for cell $n$ in batch $s$. Denote:
- $\mathcal{T}_s$ = set of proteins measured in batch $s$
- $\mathcal{I} = \bigcap_s \mathcal{T}_s$ = shared proteins across batches
- $\mathcal{U} = \bigcup_s \mathcal{T}_s$ = union of all proteins

The binary mask $M_n^{(s)}$ indicates which features were actually observed:
$$
M_{np}^{(s)} =
\begin{cases}
1 & \text{if } p \in \mathcal{T}_s \\
0 & \text{otherwise}
\end{cases}
$$

This binary mask is generated automatically when using CytoVI's 'merge_batches' function:

```
>>> adata_list = [adata_batch1, adata_batch2]
>>> adata = scvi.external.CytoVI.merge_batches(adata_list)
```

<!-- update API -->

Consecutively, the encoder uses only the shared features $\mathcal{I}$:

$$
z_n \sim q(z_n \mid x_n^{(\mathcal{I})})
$$

The decoder reconstructs $\mathcal{U}$, yielding a reconstructed protein expression vector $x_n^{\mathcal{(U)}}$ for each cell.

## Tasks
Here we provide an overview of some of the tasks that scVI can perform. Please see {class}`scvi.external.CYTOVI` for the full API reference.

### Dimensionality reduction

For dimensionality reduction, the mean of the approximate posterior $q_\eta(z_n \mid x_n, s_n)$ is returned by default.
This is achieved using the method:

```
>>> latent = model.get_latent_representation()
>>> adata.obsm["X_CytoVI"] = latent
```

Users may also return samples from this distribution, as opposed to the mean by passing the argument `give_mean=False`.
The latent representation can be used to create a nearest neighbor graph with scanpy with:

```
>>> import scanpy as sc
>>> sc.pp.neighbors(adata, use_rep="X_CytoVI")
>>> adata.obsp["distances"]
```

### Transfer learning

A CytoVI model can be pre-trained on reference data and updated with query data using {func}`~scvi.external.CYTOVI.load_query_data`, which then facilitates transfer of metadata like cell type annotations using {func}`~scvi.external.CYTOVI.impute_categories_from_reference`.

```
>>> model_query = scvi.external.CYTOVI.load_query_data(adata = adata_query, reference_model=model)
>>> model_query.is_trained = True
>>> adata_query.obs['imputed_label'] = model_query.impute_categories_from_reference(adata_reference, cat_key = 'cell_type')
```
See the {doc}`/user_guide/background/transfer_learning` guide for more information.

### Label-free differential abundance
CytoVI supports label-free differential abundance estimation via {func}`~scvi.external.CYTOVI.differential_abundance` to detect shifts in cellular composition across sample-level covariates (e.g., disease vs. control), as described in [^ref3].

This method:
- Aggregates the approximate posterior distributions $q(z \mid x_n)$ across all cells within each sample to obtain a sample-level latent density.
- Compares the average latent densities across sample groups $A$ and $B$ (e.g., condition vs. control) to compute a log-ratio score $r_{A,B}(z)$.
- Identifies enriched or depleted regions of latent space without requiring clustering.

This approach enables cluster-free detection of condition-associated cell states by directly comparing their latent representations across groups.

```
>>> da_res = model.differential_abundance(adata, groupby='group')
```

### Normalization/denoising/imputation of expression
In {func}`~scvi.model.CYTOVI.get_normalized_expression` CytoVI returns the expected value of $x_{n}^{(s)}$ under the approximate posterior. For one cell $n$, this can be written as:

$$
   \mathbb{E}_{q_\eta(z_n \mid x_n)}\left[f_x\left( z_n, s_n \right) \right]
$$

By default we decode the reconstructed protein expression across all possible batches $S$ and return the mean, yielding a batch corrected version of the protein expression.

In case of overlapping antibody panels, CytoVI by default returns protein expression of all markers $x_{n,\mathcal{U}}$, thereby effectively imputing missing proteins.

### Differential expression
Differential expression analysis is achieved with {func}`~scvi.external.CYTOVI.differential_expression`. CytoVI tests differences in magnitude of $f_x\left( z_n, s_n \right)$. More info is in {doc}`/user_guide/background/differential_expression`.

If a sample_key is provided, CytoVI by default samples equal numbers of cells for each patient for differential expression computation.

```
>>> de_res = model.differential_expression(adata, groupby='group')
```

### Data simulation
Data can be generated from the model using the posterior predictive distribution in {func}`~scvi.external.CYTOVI.posterior_predictive_sample`.
This is equivalent to feeding a cell through the model, sampling from the posterior
distributions of the latent variables, retrieving the likelihood parameters (of $p(x \mid z, s)$), and finally, sampling from this distribution.

### RNA/modality imputation
CytoVI enables cross-modal imputation by leveraging a shared latent space between datasets with overlapping protein features via {func}`~scvi.external.CYTOVI.impute_rna_from_reference`. For example, transcriptomic profiles from a CITE-seq reference can be imputed into a flow cytometry dataset.

To do this:
- A joint CytoVI model is trained on both datasets, aligning them in a shared latent space.
- For each query cell (e.g., flow cytometry), its k-nearest neighbors (default: k = 20) are identified in the latent space from the reference (e.g., CITE-seq).
- The imputed RNA expression is obtained by averaging the RNA profiles of those neighbors.

This approach allows label-free imputation of unobserved modalities, such as gene expression, in datasets where only protein measurements are available.

```
>>> adata_imputed_rna = model.impute_rna_from_reference(reference_batch='CITE_seq', adata_rna = adata_rna, layer_key='rna_normalized', return_query_only = True)
```

[^ref1]:
    Florian Ingelfinger, Nathan Levy, Can Ergen, Artemy Bakulin, Alexander Becker, Pierre Boyeau, Martin Kim, Diana Ditz, Jan Dirks, Jonas Maaskola, Tobias Wertheimer, Robert Zeiser, Corinne C. Widmer, Ido Amit, Nir Yosef,
    _CytoVI: Deep generative modeling of antibody-based single cell technologies_,
    [bioRxiv](https://doi.org/).

<!-- update ref as soon as its online -->

[^ref2]:
    Adam Gayoso\*, Zoë Steier\*, Romain Lopez, Jeffrey Regier, Kristopher L Nazor, Aaron Streets, Nir Yosef (2021),
    _Joint probabilistic modeling of single-cell multi-omic data with totalVI_,
    [Nature Methods](https://www.nature.com/articles/s41592-020-01050-x).

[^ref3]:
     Pierre Boyeau, Justin Hong, Adam Gayoso, Martin Kim, Jose L McFaline-Figueroa, Michael Jordan, Elham Azizi, Can Ergen, Nir Yosef (2024),
    _Deep generative modeling of sample-level heterogeneity in single-cell genomics_,
    [bioRxiv](https://doi.org/10.1101/2022.10.04.510898).
