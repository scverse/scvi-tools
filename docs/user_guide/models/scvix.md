# scVI-X

**scVI-X** [^ref1] (single-cell Variational Inference — cross-technology;
Python class {class}`~scvi.external.SCVIX`) is a deep generative model for
learning technology-invariant cell-state representations across diverse
single-cell sequencing platforms. It extends scVI [^ref2] with three targeted
architectural changes that together enable strong cross-technology integration
and zero-shot query mapping without fine-tuning.

The key advantages of scVI-X over scVI are:

- **Cross-technology integration**: Explicitly models the sequencing assay
  (e.g., 10x single-cell vs. single-nucleus) as a distinct source of
  variation, yielding substantially better mixing across platforms.
- **Zero-shot query mapping**: Trained models can embed new cells from
  previously unseen arrays via a single forward pass, with no retraining
  required.
- **Assay-level adversarial training**: Adversarial correction is applied at
  the assay level rather than the batch level, reducing the risk of
  overintegration when cell-type compositions differ between batches.
- **Assay-aware decoding**: Learned gene-level assay biases in the output
  layer capture technology-specific transcript-detection differences and
  provide interpretable, reproducible signatures across datasets.
- **Variational batch embeddings**: Optionally replaces one-hot batch
  encodings with low-dimensional variational embeddings, reducing parameter
  count and improving identifiability through a KL regularisation term.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/scrna/scvix_tutorial`
```

## Preliminaries

scVI-X takes as input a scRNA-seq count matrix $X$ with $N$ cells and $G$
genes. For each cell $n$ it requires:

- An **assay label** $a_n$ identifying the sequencing technology or
  suspension type (e.g., 10x Chromium single-cell, 10x Chromium
  single-nucleus, Smart-seq2). This is the primary covariate that scVI-X is
  designed to correct for.
- A **batch label** $b_n$ for experimental nuisance variation within an
  assay (e.g., sequencing run, donor, processing day).

Optionally, scVI-X can also accept:

- Cell-type labels $y_n$ for supervised integration via a cell-type
  biased prior.
- An **adversarial group label** $g_n$ (defaults to $a_n$) that defines
  which grouping the adversarial classifier acts on. To reduce
  overintegration if cell-type composition varies across assays.
- Additional categorical or continuous covariates.

## Generative process

scVI-X posits a conditional VAE generative model. For each cell $n$:

1. A latent cell-state vector $z_n$ is drawn from the prior:

    $z_n \sim p(z)$

    where $p(z)$ is one of: a standard Gaussian, a mixture of Gaussians
    (MOG), a variational amortized mixture of posteriors (VampPrior).

2. The normalised expression scale $h_n$ is decoded from $z_n$ conditioned
   on the batch $b_n$ and assay $a_n$:

    $h_n = \mathrm{softmax}\!\left(f_\theta(z_n,\, b_n) + \delta_{a_n}\right)$

    where $f_\theta$ is a fully-connected decoder network, and
    $\delta_{a_n} \in \mathbb{R}^G$ is a learned, assay-specific output bias
    vector that captures gene-level differences in transcript-detection
    efficiency across sequencing platforms.

3. Observed counts are generated from a negative binomial (or ZINB / Poisson)
   likelihood:

    $x_{ng} \mid h_{ng} \sim \mathrm{NegativeBinomial}(l_n\, h_{ng},\; r_{ng})$

    where $l_n$ is the observed library size and $r_{ng}$ is the
    gene-specific (or cell-specific) inverse dispersion.

The latent variables are summarised below:

```{eval-rst}
.. list-table::
   :widths: 20 80 15
   :header-rows: 1

   * - Variable
     - Description
     - Code name
   * - :math:`z_n \in \mathbb{R}^L`
     - Technology-invariant cell-state embedding.
     - ``z``
   * - :math:`h_n \in \mathbb{R}^G`
     - Normalised gene expression scale.
     - ``px_scale``
   * - :math:`l_n \in \mathbb{R}^+`
     - Observed library size (log-transformed).
     - ``library``
   * - :math:`r_{ng} \in \mathbb{R}^+`
     - Gene- (and optionally cell-) specific inverse dispersion.
     - ``px_r``
   * - :math:`\delta_{a_n} \in \mathbb{R}^G`
     - Assay-specific gene-level output bias.
     - output layer weights conditioned on ``assay_index``
```

### Variational batch embeddings (optional)

When `batch_representation="embedding"` is selected, each batch $b$ is
represented by a low-dimensional embedding
$e_b \sim \mathcal{N}(\mu_b, \sigma_b^2 I)$ rather than a one-hot vector.
This reduces the number of decoder parameters. This embedding can be variational
and then adds a KL regularisation
term to the loss that improves identifiability of batch representations:

$\mathcal{L}_\text{embed} = \mathrm{KL}\!\left(q(e_b) \,\|\, \mathcal{N}(0,I)\right)$

## Inference

scVI-X uses amortised variational inference. The approximate posterior is:

$q_\phi(z_n \mid x_n, a_n) = \mathcal{N}\!\left(\mu_\phi(x_n, a_n),\; \sigma^2_\phi(x_n, a_n)\, I\right)$

where $\mu_\phi$ and $\sigma^2_\phi$ are encoder neural networks. The assay
index $a_n$ is injected into the encoder through **conditional layer
normalisation** (controlled by `conditional_norm=True`, the default): instead
of a single global scale and shift applied after each hidden layer, the model
learns per-assay scale and shift parameters $\gamma_{a_n}$ and $\beta_{a_n}$.
This provides a lightweight but effective mechanism for correcting
assay-specific distributional shifts in the encoder.

The evidence lower bound (ELBO) trained is:

$\mathcal{L} = \mathbb{E}_{q_\phi(z_n \mid x_n, a_n)}\!\left[\log p_\theta(x_n \mid z_n, b_n, a_n)\right] - \mathrm{KL}\!\left(q_\phi(z_n \mid x_n, a_n) \,\|\, p(z_n)\right) + \mathcal{L}_\text{embed}$

### Adversarial training

When multiple assays are present, scVI-X adds an adversarial
classifier that acts on the latent space $z_n$. The classifier is trained to
predict the adversarial label $a_n$ (typically the assay), while the
encoder is simultaneously trained to fool it:

$\mathcal{L}_\text{adv} = -\mathbb{E}\left[\log p_\psi(\hat{a}_n \neq a_n \mid z_n)\right]$

By conditioning adversarial correction at the **assay** level rather than the
fine-grained batch level, scVI-X reduces the risk of overintegration in
settings where cell-type compositions differ between donors or conditions
within the same assay.

Adversarial training is automatically enabled when more than one assay is
present; the adversarial label can be set to any registered categorical
covariate via `adversarial_key`.

## Prior options

scVI-X supports four prior distributions for $z_n$, selected via the `prior`
argument at model initialisation. We found Gaussian to perform well in all tested scenarios.

- **`'gaussian'`** (default): Standard normal $\mathcal{N}(0, I)$. Fastest
  to train; suitable for most use cases.
- **`'mog'`**: Mixture of Gaussians with $K$ learnable components. Encourages
  a more structured latent space.
- **`'vamp'`**: VampPrior ([Tomczak & Welling, 2018](https://doi.org/10.48550/arXiv.1705.07120)).
  Prior modes are anchored to learned pseudoinputs drawn from the data,
  making the prior data-adaptive.
- **`'mog_celltype'`**: Mixture of Gaussians with one component per cell-type
  label. Requires `labels_key` to be set in `setup_anndata`. Guides
  integration by biasing the prior toward annotated cell-type clusters.

## Tasks

### Latent representation

The primary output of scVI-X is a low-dimensional, technology-corrected
embedding of each cell:

```python
import scvi

scvi.external.SCVIX.setup_anndata(
    adata,
    batch_key="batch",
    assay_key="assay",  # sequencing technology / suspension type
)
model = scvi.external.SCVIX(
    adata
)  # defaults: embedding, n_latent=20, n_layers=2, n_hidden=512
model.train()  # defaults: batch_size=1024, n_epochs_kl_warmup=5
adata.obsm["X_scVIX"] = model.get_latent_representation()
```

### Normalised expression

Batch- and assay-corrected normalised expression can be obtained by decoding
from the latent space while fixing the batch and/or assay to a reference
value. This is analogous to `get_normalized_expression` in scVI and supports
the `transform_batch` and `transform_assay` arguments:

```python
# Expression as it would look in a given reference batch
norm_expr = model.get_normalized_expression(transform_batch="batch_1")
```

### Zero-shot query mapping

scVI-X supports query-to-reference mapping via the scArches framework without
any retraining. A query dataset with new batches (but the same gene set and assay) can
be embedded directly without additional training.

```python
# Save the reference model
model.save("scvix_reference/")

# Prepare the query and embed zero-shot
scvi.external.SCVIX.prepare_query_anndata(adata_query, "scvix_reference/")
query_model = scvi.external.SCVIX.load_query_data(adata_query, "scvix_reference/")
query_model.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})
adata_query.obsm["X_scVIX"] = query_model.get_latent_representation()
```

### Differential expression

Standard differential expression analysis between cell groups is available
via the inherited `differential_expression` method:

```python
de_results = model.differential_expression(groupby="cell_type", group1="T cell")
```

### Batch embeddings

When `batch_representation="embedding"` is used, the learned batch embeddings
can be retrieved and used for downstream analysis (e.g., visualising
batch-level variation or clustering batches):

```python
model = scvi.external.SCVIX(adata, batch_representation="embedding")
model.train()
batch_emb = model.get_batch_representation()  # shape: (n_cells, embedding_dim)
```

## Key setup_anndata arguments

```{eval-rst}
.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Argument
     - Description
   * - ``batch_key``
     - Column in ``adata.obs`` for experimental batch (donor, run, etc.).
   * - ``assay_key``
     - Column in ``adata.obs`` for sequencing technology / suspension type.
       Drives conditional layer normalisation and assay-output biases.
```

[^ref1]:
    Can Ergen, Ori Kronfeld, Martin Kim, Shiyi Yang, Florian Ingelfinger,
    Nir Yosef (in preparation),
    _scvi-X: Learning technology-invariant cell states_.

[^ref2]:
    Romain Lopez, Jeffrey Regier, Michael B Cole, Michael I Jordan, Nir Yosef
    (2018),
    _Deep generative modeling for single-cell transcriptomics_,
    [Nature Methods](https://doi.org/10.1038/s41592-018-0229-2).
