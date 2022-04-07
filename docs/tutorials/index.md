# Tutorials

The easiest way to get familiar with scvi-tools is to follow along with our tutorials.
Many are also designed to work seamlessly in Google Colab, a free cloud computing platform.
Tutorials by default work with the latest installable version of scvi-tools. To view older tutorials,
change the documentation version using the tab at the bottom of the left sidebar.

:::{note}
For questions about using scvi-tools, or broader questions about modeling data, please use our [forum]. Checkout the [ecosystem] for additional models powered by scvi-tools.
:::

## Quick start

```{eval-rst}
.. nbgallery::

   notebooks/api_overview
   notebooks/data_loading
   notebooks/python_in_R

```

## User tutorials

### scRNA-seq

```{eval-rst}
.. nbgallery::

   notebooks/harmonization
   notebooks/scvi_in_R
   notebooks/tabula_muris
   notebooks/scarches_scvi_tools
   notebooks/query_hlca_knn
   notebooks/seed_labeling
   notebooks/linear_decoder
   notebooks/AutoZI_tutorial
   notebooks/cellassign_tutorial
   notebooks/amortized_lda

```

### ATAC-seq

```{eval-rst}
.. nbgallery::

   notebooks/PeakVI
   notebooks/peakvi_in_R

```

### Spatial transcriptomics

```{eval-rst}
.. nbgallery::

   notebooks/DestVI_tutorial
   notebooks/DestVI_in_R
   notebooks/gimvi_tutorial

```

### Multimodal

```{eval-rst}
.. nbgallery::

   notebooks/totalVI
   notebooks/cite_scrna_integration_w_totalVI
   notebooks/totalVI_reference_mapping
   notebooks/totalvi_in_R
   notebooks/MultiVI_tutorial

```

## Contributed tutorials

```{eval-rst}
.. nbgallery::

   notebooks/stereoscope_heart_LV_tutorial
   notebooks/cell2location_lymph_node_spatial_tutorial

```

## Developer tutorials

Here we feature tutorials to guide you through the construction of a model with scvi-tools. For an example of how scvi-tools can be used in an independent package, see our GitHub [template].

:::{note}
For questions about developing with scvi-tools, please use our [forum] or [gitter].
:::

```{toctree}
:maxdepth: 1

notebooks/data_tutorial
notebooks/module_user_guide
notebooks/model_user_guide
```

### Glossary

```{eval-rst}
.. tab-set::

    .. tab-item:: Model

        A Model class inherits :class:`~scvi.model.base.BaseModelClass` and is the user-facing object for interacting with a module.
        The model has a `train` method that learns the parameters of the module, and also contains methods
        for users to retrieve information from the module, like the latent representation of cells in a VAE.
        Conventionally, the post-inference model methods should not store data into the AnnData object, but
        instead return "standard" Python objects, like numpy arrays or pandas dataframes.


    .. tab-item:: Module

        A module is the lower-level object that defines a generative model and inference scheme. A module will
        either inherit :class:`~scvi.module.base.BaseModuleClass` or :class:`~scvi.module.base.PyroBaseModuleClass`.
        Consequently, a module can either be implemented with PyTorch alone, or Pyro. In the PyTorch only case, the
        generative process and inference scheme are implemented respectively in the `generative` and `inference` methods,
        while the `loss` method computes the loss, e.g, ELBO in the case of variational inference.
```

```{eval-rst}
.. tab-set::

    .. tab-item:: TrainingPlan

        The training plan is a PyTorch Lightning Module that is initialized with a scvi-tools module object.
        It configures the optimizers, defines the training step and validation step, and computes metrics to be
        recorded during training. The training step and validation step are functions that take data, run it through
        the model and return the loss, which will then be used to optimize the model parameters in the Trainer.
        Overall, custom training plans can be used to develop complex inference schemes on top of modules.



    .. tab-item:: Trainer

        The :class:`~scvi.train.Trainer` is a lightweight wrapper of the PyTorch Lightning Trainer. It takes as input
        the training plan, a training data loader, and a validation dataloader. It performs the actual training loop, in
        which parameters are optimized, as well as the validation loop to monitor metrics. It automatically handles moving
        data to the correct device (CPU/GPU).
```

[ecosystem]: https://scvi-tools.org/ecosystem
[forum]: https://discourse.scvi-tools.org/
[gitter]: https://gitter.im/scvi-tools/development
[template]: https://github.com/YosefLab/scvi-tools-skeleton
