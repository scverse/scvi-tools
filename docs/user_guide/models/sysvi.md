# SysVI

**sysVI** (cross-SYStem Variational Inference,
Python class {class}`~scvi.external.SysVI`)
is a representation learning models that can remove substantial batch effects.

The advantages of SysVI are:

-   Improved integration: For datasets with **substantial batch effects**
(e.g., cross-species or organoid-tissue), where other models often fail.
It provides a good tradeoff between batch correction and preservation of
cell-type and sub-cell-type biological variation.
- Tunable integration: The **integration strength is directly tunable**
via cycle consistency loss.
- Generally applicable: The model operates on
**approximately normally distributed data**
(e.g. normalized and log+1 transformed scRNA-seq data), which makes it
more generally applicable than just scRNA-seq.
- Scalable: Can integrate very large datasets if using a GPU.

The limitations of SysVI include:

-   Weak batch effects: For datasets with **small batch effects**
(e.g. multiple subjects from a single laboratory) we recommend using scVI instead,
as it has slightly higher biological preservation in this setting.
For determining whether a dataset has substantial batch effects
please refer to our paper.
- Model selection: The best performance is achieved if
**selecting the best model** from multiple
runs with a few different cycle consistency loss weights and random seed
initialisations, as explained in the tutorial.
However, we provide **defaults** that generate decent results in
many settings.


```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/scrna/sysVI`
```

```{topic} References:

-  Paper: Hrovatin and Moinfar, et al.
Integrating single-cell RNA-seq datasets with substantial batch effects.
bioRxiv (2023): https://doi.org/10.1101/2023.11.03.565463
- Talk on caveats of scRNA-seq integration and strategies for removing
substantial batch effects: https://www.youtube.com/watch?v=i-a4BjAn90E
```

## Method background

The model is based on a variational autoencoder (VAE), with the integrated
representation corresponding to the latent space embedding of the cells.

### Stronger batch correction with cycle-consistency loss

Vanilla VAEs struggle to achieve strong batch correction without loosing
substantial biological variation. This issue arises as the VAE loss
does not directly penalize the presence of batch covariate information in the
latent space.
Instead, conditional VAEs assume that batch covariate information will be
omitted from the latent space, which has limited-capacity,
as it is separately injected into the decoder. Namely, its presence in the
latent space is "unnecessary" for the reconstruction (Hrovatin and Moinfar, 2023).

To achieve stronger integration than vanilla VAEs, SysVI employs
cycle-consistency loss in the latent space. In particular, the model embeds a cell
from one system (i.e. the covariate representing substantial batch effect)
into latent space and then decodes it using another category of the system covariate.
In this way it generates a biologically identical cell with a
different batch effect. The generated cell is then likewise embedded into the
latent space and the distance between the embeddings of the original and
the switched-batch cell are computed. The model is trained to minimize this distance.

:::{figure} figures/sysvi_cycleconsistency.png
:align: center
:alt: Cycle consistency loss used to increase batch correction in SysVI.
:class: img-fluid
:::

Benefits of this approach:
- As only cells with identical biological background are compared, this method
retains good biological preservation even when removing
substantial batch effects. This distinguishes it from alternative approaches
that compare cells with different biological backgrounds
(e.g. via adversarial loss; see Hrovatin and Moinfar (2023) for details).
- The integration strength can be directly tuned via the cycle-consistency
loss weight.

### Improved biological preservation via the VampPrior

Vanilla VAEs employ standard normal prior for regularizing latent space.
However, this prior is very restrictive and can lead to loss of
important biological variation in the latent space.

Instead, we use the
VampPrior ([Tomczak, 2017](https://doi.org/10.48550/arXiv.1705.07120)),
which permits a more expressive latent space. VampPrior is a multi-modal
prior for which the mode positions are learned during the training.

:::{figure} figures/sysvi_vampprior.png
:align: center
:alt: VampPrior used to increase the preservation of biological variation in SysVI.
:class: img-fluid
:::

Benefits of this approach:
- More expressive latent space leads to increased preservation of
biological variability.
- VampPrior was more robust with respect to the number of modes than the
better-know Gaussian mixture prior.

### Application flexibility due to using normally distributed inputs

Many scRNA-seq integration models are specially designed to work with
scRNA-seq data, e.g. raw counts that follows negative binomial distribution.
However, due to this, these models can not be directly used for other
types of data.

We observed that for representation learning this specialised setup is not
strictly required. - SysVI is designed for data following normal distribution, while
performing competitively in comparison to the more specialised models
on scRNA-seq data.
To make scRNA-seq data approximately normally distributed we preprocess it via
size-factor normalization and log+1 transformation.

Thus, SysVI could be also applied to other types of normally distributed data.
However, we did not specifically test its performance on other data types.

## Other tips & tricks for data integration

Besides the benefits of the SysVI model, our paper
([Hrovatin and Moinfar, 2023](https://doi.org/10.1101/2023.11.03.565463))
and
[talk](https://www.youtube.com/watch?v=i-a4BjAn90E)
provide additional advice on scRNA-seq integration that apply beyond SysVI.
The two most important insights are:
- Try to make the **integration task as easy for the model** as possible.
This means that data should be pre-processed in a way that already eliminates
some of the batch differences, when possible:
  - Use intersection of HVGs across batches with substantial batch effects
  (e.g. the systems).
  - Mitigate known technical artefacts, such as ambient gene expression
  ([Hrovatin and Sikkema, 2024](https://doi.org/10.1038/s41592-024-02532-y)).
- Ensure that **the metrics used to evaluate integration are of high-quality**:
  - They should be able to capture the key properties required for downstream tasks.
  For example, the standard cell-type based biological preservation metrics do
  not assess whether subtler biological differences, such as within-cell-type
  disease effects, are preserved.
  - Be cautious of potential biases within integration metric scores. -
  The scores may not directly correspond to the desired data property,
  being influenced by other factors, or
  certain models may be able to trick the metrics.
