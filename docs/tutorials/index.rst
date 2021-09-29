Tutorials
==========

The easiest way to get familiar with scvi-tools is to follow along with our tutorials.
Many are also designed to work seamlessly in Google Colab, a free cloud computing platform.
Tutorials are designed to work with the latest installable version of scvi-tools.

.. note:: For questions about using scvi-tools, or broader questions about modeling data, please use our forum_.

.. _forum: https://discourse.scvi-tools.org/


Quick start
-----------

.. nbgallery::

   notebooks/api_overview
   notebooks/data_loading
   notebooks/scvi_in_R


User tutorials
--------------

scRNA-seq
^^^^^^^^^^

.. nbgallery::

   notebooks/harmonization
   notebooks/seed_labeling
   notebooks/linear_decoder
   notebooks/scarches_scvi_tools
   notebooks/AutoZI_tutorial
   notebooks/cellassign_tutorial


ATAC-seq
^^^^^^^^

.. nbgallery::

   notebooks/PeakVI
   notebooks/peakvi_in_R


Spatial transcriptomics
^^^^^^^^^^^^^^^^^^^^^^^


.. nbgallery::

   notebooks/DestVI_tutorial
   notebooks/gimvi_tutorial


Multimodal
^^^^^^^^^^^

.. nbgallery::

   notebooks/totalVI
   notebooks/cite_scrna_integration_w_totalVI
   notebooks/totalVI_reference_mapping
   notebooks/totalvi_in_R
   notebooks/MultiVI_tutorial



Contributed tutorials
---------------------

.. nbgallery::

   notebooks/scVI_DE_worm
   notebooks/stereoscope_heart_LV_tutorial


Developer tutorials
---------------------

Here we feature tutorials to guide you through the construction of a model with scvi-tools. For an example of how scvi-tools can be used in an independent package, see our GitHub template_.

.. note:: For questions about developing with scvi-tools, please use our forum_ or gitter_.

.. _forum: https://discourse.scvi-tools.org/
.. _gitter: https://gitter.im/scvi-tools/development
.. _template: https://github.com/YosefLab/scvi-tools-skeleton


.. toctree::
   :maxdepth: 1

   notebooks/data_tutorial
   notebooks/module_user_guide
   notebooks/model_user_guide


Glossary
^^^^^^^^^

.. tabs::

    .. tab:: Model

        A Model class inherits :class:`~scvi.model.base.BaseModelClass` and is the user-facing object for interacting with a module.
        The model has a `train` method that learns the parameters of the module, and also contains methods
        for users to retrieve information from the module, like the latent representation of cells in a VAE.
        Conventionally, the post-inference model methods should not store data into the AnnData object, but
        instead return "standard" Python objects, like numpy arrays or pandas dataframes.

    .. tab:: Module

        A module is the lower-level object that defines a generative model and inference scheme. A module will
        either inherit :class:`~scvi.module.base.BaseModuleClass` or :class:`~scvi.module.base.PyroBaseModuleClass`.
        Consequently, a module can either be implemented with PyTorch alone, or Pyro. In the PyTorch only case, the
        generative process and inference scheme are implemented respectively in the `generative` and `inference` methods,
        while the `loss` method computes the loss, e.g, ELBO in the case of variational inference.

.. tabs::

    .. tab:: TrainingPlan


        The training plan is a PyTorch Lightning Module that is initialized with a scvi-tools module object.
        It configures the optimizers, defines the training step and validation step, and computes metrics to be
        recorded during training. The training step and validation step are functions that take data, run it through
        the model and return the loss, which will then be used to optimize the model parameters in the Trainer.
        Overall, custom training plans can be used to develop complex inference schemes on top of modules.


    .. tab:: Trainer

        The :class:`~scvi.train.Trainer` is a lightweight wrapper of the PyTorch Lightning Trainer. It takes as input
        the training plan, a training data loader, and a validation dataloader. It performs the actual training loop, in
        which parameters are optimized, as well as the validation loop to monitor metrics. It automatically handles moving
        data to the correct device (CPU/GPU).
