# CytoVI

**CytoVI** (Python class {class}`~scvi.external.CYTOVI`) is a generative model for cytometry data that leverages deep probabilistic latent variable modeling to enable denoising, imputation, integration, and differential analysis across technologies and batches.

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
- {doc}`/tutorials/notebooks/cytometry/CytoVI_tutorial` # note: to populate
```

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
To integrate datasets with different protein panels, CytoVI employs a masking strategy inspired by {class}`~scvi.model.TOTALVI`.

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

### Transfer learning

### Label-free differential abundance

### Normalization/denoising/imputation of expression

### Differential expression

### Data simulation

### RNA/modality imputation


