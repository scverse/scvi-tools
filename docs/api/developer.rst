=========
Developer
=========


Import scvi-tools as::

   import scvi


.. currentmodule:: scvi


Data Loaders
------------

.. currentmodule:: scvi

DataLoaders for loading tensors from AnnData objects. DataSplitters for splitting data into train/test/val.


.. autosummary::
   :toctree: reference/
   :nosignatures:

   dataloaders.AnnDataLoader
   dataloaders.AnnTorchDataset
   dataloaders.ConcatDataLoader
   dataloaders.SemiSupervisedDataLoader
   dataloaders.DataSplitter
   dataloaders.SemiSupervisedDataSplitter


Distributions
-------------

.. currentmodule:: scvi

Parameterizable probability distributions.


.. autosummary::
   :toctree: reference/
   :nosignatures:

   distributions.NegativeBinomial
   distributions.NegativeBinomialMixture
   distributions.ZeroInflatedNegativeBinomial


Model (Base)
------------

.. currentmodule:: scvi

These classes should be used to construct user-facing model classes.

.. autosummary::
    :toctree: reference/
    :nosignatures:

    model.base.BaseModelClass
    model.base.VAEMixin
    model.base.RNASeqMixin
    model.base.ArchesMixin
    model.base.UnsupervisedTrainingMixin

Module
------

.. currentmodule:: scvi

Existing module classes with respective generative and inference procedures.

.. autosummary::
   :toctree: reference/
   :nosignatures:

   module.VAE
   module.LDVAE
   module.TOTALVAE
   module.SCANVAE
   module.JVAE
   module.AutoZIVAE
   module.Classifier

Module (Base)
-------------

.. currentmodule:: scvi

These classes should be used to construct module classes that define generative models and inference schemes.

.. autosummary::
   :toctree: reference/
   :nosignatures:

   module.base.LossRecorder
   module.base.BaseModuleClass
   module.base.PyroBaseModuleClass
   module.base.auto_move_data
   

Neural networks
---------------

.. currentmodule:: scvi

Basic neural network building blocks.

.. autosummary::
   :toctree: reference/
   :nosignatures:
   
   nn.FCLayers
   nn.Encoder
   nn.Decoder
   nn.one_hot

Train
-----

.. currentmodule:: scvi


TrainingPlans define train/test/val optimization steps for modules.

.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst
   :nosignatures:

   train.TrainingPlan
   train.SemiSupervisedTrainingPlan
   train.AdversarialTrainingPlan
   train.PyroTrainingPlan
   train.Trainer
   train.TrainRunner

Utilities
---------

.. currentmodule:: scvi

Utility functions used by scvi-tools.

.. autosummary::
   :toctree: reference/
   :nosignatures:

   utils.DifferentialComputation
   utils.track
