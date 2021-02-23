=========
Developer
=========


Import scvi-tools as::

   import scvi


.. currentmodule:: scvi

Compose
-------

.. currentmodule:: scvi


Architectures
~~~~~~~~~~~~~

.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst

   compose.FCLayers
   compose.Encoder
   compose.Decoder

Module classes
~~~~~~~~~~~~~~

.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst

   compose.LossRecorder
   compose.BaseModuleClass
   compose.PyroBaseModuleClass

Functions
~~~~~~~~~

.. autosummary::
   :toctree: reference/

   compose.one_hot
   compose.auto_move_data

Distributions
-------------

.. currentmodule:: scvi


.. autosummary::
   :toctree: reference/

   distributions.NegativeBinomial
   distributions.NegativeBinomialMixture
   distributions.ZeroInflatedNegativeBinomial


Model (Base)
------------

.. currentmodule:: scvi


.. autosummary::
    :toctree: reference/

    model.base.BaseModelClass
    model.base.VAEMixin
    model.base.RNASeqMixin
    model.base.ArchesMixin
    model.base.TrainRunner

Modules
-------

.. currentmodule:: scvi


.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst

   modules.VAE
   modules.LDVAE
   modules.TOTALVAE
   modules.SCANVAE
   modules.JVAE
   modules.AutoZIVAE
   modules.Classifier

Data Loaders
------------

.. currentmodule:: scvi


.. autosummary::
   :toctree: reference/

   dataloaders.AnnDataLoader
   dataloaders.AnnTorchDataset
   dataloaders.ConcatDataLoader
   dataloaders.SemiSupervisedDataLoader
   dataloaders.DataSplitter
   dataloaders.SemiSupervisedDataSplitter

Lightning
---------

.. currentmodule:: scvi


PyTorch lightning is used to train our modules. TrainingPlans define train/test/val optimization
steps for modules like `TOTALVAE`, `SCANVAE`, etc.

.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst

   lightning.TrainingPlan
   lightning.SemiSupervisedTrainingPlan
   lightning.AdversarialTrainingPlan
   lightning.PyroTrainingPlan
   lightning.Trainer

Utilities
---------

.. currentmodule:: scvi


.. autosummary::
   :toctree: reference/

   utils.DifferentialComputation
   utils.track
