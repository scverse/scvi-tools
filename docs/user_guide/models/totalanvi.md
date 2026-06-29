# TotalANVI

**TotalANVI** {cite:p}`GayosoSteier21` (Python class {class}`~scvi.external.TOTALANVI`) is a semi-supervised generative model of CITE-seq RNA and protein data.
Similar to how {doc}`scANVI <scanvi>` extends {doc}`scVI <scvi>`, TotalANVI extends {doc}`TotalVI <totalvi>` by leveraging partial cell type annotations
to jointly model gene expression and protein abundance, predict labels for unannotated cells, and impute missing protein measurements.

The advantages of TotalANVI are:

-   Unified multimodal (RNA + protein) semi-supervised framework.
-   Supports training from scratch or fine-tuning a pretrained TotalVI model.
-   Handles cells with missing protein panels.
-   Scalable to very large datasets (>1 million cells).

The limitations of TotalANVI include:

-   Effectively requires a GPU for fast inference.
-   Requires at least a subset of cells to have cell type annotations.
-   May not scale to a very large number of cell types.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/multimodal/totalanvi`

```

## Preliminaries

Let $X \in \mathbb{N}^{N \times G}$ be the RNA count matrix for $N$ cells and $G$ genes, and $Y \in \mathbb{N}^{N \times T}$ be the protein count matrix for $T$ proteins. Not all cells are required to have measured proteins. Let $S \in \mathbb{R}^{N \times P}$ denote observed covariates (e.g., batch). Let $c_n \in \{1, \ldots, K, \text{unlabeled}\}$ denote the (partial) cell type annotation for cell $n$, where cells assigned the `unlabeled_category` have unknown cell type.

## Generative process

TotalANVI combines TotalVI's protein–RNA generative model with scANVI's cell-type-conditioned hierarchical latent space.

For each cell $n$, the generative process is:

$$
c_n \sim \text{Cat}(\pi)
$$

$$
u_n \mid c_n \sim \mathcal{N}(0, I)
$$

$$
z_n \mid u_n, c_n \sim \mathcal{N}\!\left(\mu_\phi(u_n, c_n),\, \sigma^2_\phi(u_n, c_n) \cdot I\right)
$$

$$
\ell_n \sim \text{LogNormal}(\ell_\mu, \ell_{\sigma^2})
$$

$$
\rho_n \propto f_w(z_n, s_n), \quad x_n \mid z_n, \ell_n, s_n \sim \text{NB}\!\left(\ell_n \rho_n,\, \theta\right)
$$

$$
\mu_{y_n} = f_h(z_n, s_n)
$$

$$
v_{n,j} \mid z_n \sim \text{Bernoulli}(\pi_{n,j}^{\text{fore}})
$$

$$
y_{n,j} \mid z_n, v_{n,j}, s_n \sim \begin{cases}
\text{NB}\!\left(\alpha_{n,j}^{\text{back}},\, \beta_{n,j}\right) & \text{if } v_{n,j} = 0 \\
\text{NB}\!\left(\mu_{y_{n,j}} + \alpha_{n,j}^{\text{back}},\, \beta_{n,j}^{\text{fore}}\right) & \text{if } v_{n,j} = 1
\end{cases}
$$

where $\pi$ is a learned (or uniform) prior over cell types, $u_n$ is a cell-type-specific latent variable, $z_n$ is the cell state conditioned on $u_n$ and $c_n$, $\ell_n$ is the RNA library size, and $v_{n,j}$ indicates protein foreground. The protein background parameters $\alpha^{\text{back}}$ and $\beta$ are shared with TotalVI.

The latent variables are summarized below.

| Latent variable | Code variable | Description |
| :--- | :--- | :--- |
| $c$ | `predicted_label` | Cell type label |
| $z_1$ (= $z$) | `z` | Observation-level cell state |
| $u$ (= $z_2$) | `untran_z` | Cell-type-specific latent state |
| $\ell$ | `library` | RNA library size |
| $\beta$ | `back_rate` | Protein background rate |

## Inference

TotalANVI uses a structured variational posterior:

$$
q(z_n, u_n, \ell_n, c_n \mid x_n, y_n) = q(z_n \mid x_n, y_n)\, q(u_n \mid z_n, c_n)\, q(\ell_n \mid x_n)\, q(c_n \mid z_n)
$$

where:

- $q(z_n \mid x_n, y_n)$ is the TotalVI encoder that takes concatenated RNA and protein inputs.
- $q(u_n \mid z_n, c_n)$ is an auxiliary encoder conditioned on the current cell state and label.
- $q(\ell_n \mid x_n)$ is the RNA library size encoder.
- $q(c_n \mid z_n)$ is the classifier network that predicts cell type probabilities from the latent state.

## Training details

TotalANVI is trained by maximizing a modified ELBO that handles labeled and unlabeled cells differently.

For **unlabeled** cells, the cell type $c_n$ is marginalized out:

$$
\mathcal{L}_U = \mathbb{E}_{q(z, u, \ell)}\!\left[\sum_{c} q(c \mid z)\, \log p(x, y, z, u, \ell \mid c)\right] - \mathrm{KL}\!\left(q(c \mid z) \,\|\, p(c)\right)
$$

For **labeled** cells, the known label is used directly, and a cross-entropy classification term is added:

$$
\mathcal{L}_S = \mathcal{L}_U\big|_{c = c^*} + \alpha \cdot \mathrm{CE}\!\left(q(c \mid z),\, c^*\right)
$$

where $\alpha$ is the `classification_ratio` hyperparameter (default 50; tutorial uses 200 for large datasets).

TotalANVI can be trained from scratch via {meth}`~scvi.external.TOTALANVI.setup_anndata` / {meth}`~scvi.external.TOTALANVI.setup_mudata`, or fine-tuned from a pretrained TotalVI model via {meth}`~scvi.external.TOTALANVI.from_totalvi_model`. In the fine-tuning workflow, the protein background parameters are frozen and only the cell-type-related components are updated.

## Tasks

TotalANVI supports all tasks available in {doc}`TotalVI <totalvi>` (dimensionality reduction, normalization, protein denoising, differential expression, data simulation), plus:

### Cell type label prediction

Unannotated cells are assigned labels by taking the argmax of the classifier posterior $q(c \mid z)$:

```python
predictions = model.predict()
soft_predictions = model.predict(soft=True)  # returns DataFrame of probabilities
```

Integrated gradients attributions for interpretability are also available:

```python
predictions, attributions = model.predict(ig_interpretability=True)
```
