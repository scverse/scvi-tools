# LDVAE

**LDVAE** [^ref1] (Linearly-decoded Variational Auto-encoder, also called Linear scVI; Python class {class}`~scvi.model.LinearSCVI`)
is a flavor of scVI with a linear decoder.

The advantages of LDVAE are:

-   Can be used to interpret latent dimensions with factor loading matrix.
-   Scalable to very large datasets (>1 million cells).

The limitations of LDVAE include:

-   Less capacity than scVI, which uses a neural network decoder.
-   Less capable of integrating data with complex batch effects.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/scrna/linear_decoder`
```

## Contrasting with scVI

Here we discuss the differences between LDVAE and scVI.

-   In LDVAE, $f_w(z_n, s_n)$ is a linear function, and thus can be represented by a matrix $W$ of dimensions $G$ (genes) by $(d + k)$ (latent space dim plus covariate categories).
-   This matrix $W$ can be accessed using {func}`~scvi.model.LinearSCVI.get_loadings`
-   LDVAE does not offer transfer learning capabilities currently.

[^ref1]:
    Valentine Svensson, Adam Gayoso, Nir Yosef, Lior Pachter (2020),
    _Interpretable factor models of single-cell RNA-seq via variational autoencoders_,
    [Bioinformatics](https://academic.oup.com/bioinformatics/article/36/11/3418/5807606).
