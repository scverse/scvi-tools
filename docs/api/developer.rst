=========
Developer
=========


Import scvi-tools as::

   import scvi


.. currentmodule:: scvi


Data Loaders
------------

.. currentmodule:: scvi


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

Modules
-------

.. currentmodule:: scvi

.. autosummary::
   :toctree: reference/
   :nosignatures:

   modules.VAE
   modules.LDVAE
   modules.TOTALVAE
   modules.SCANVAE
   modules.JVAE
   modules.AutoZIVAE
   modules.Classifier

Modules (Base)
--------------

.. currentmodule:: scvi

These classes should be used to construct module classes that define generative models and inference schemes.

.. autosummary::
   :toctree: reference/
   :nosignatures:

   modules.base.LossRecorder
   modules.base.BaseModuleClass
   modules.base.PyroBaseModuleClass
   modules.base.auto_move_data
   

Neural networks
---------------

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


TrainingPlans define train/test/val optimization
steps for modules like `TOTALVAE`, `SCANVAE`, etc.

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


.. autosummary::
   :toctree: reference/
   :nosignatures:

   utils.DifferentialComputation
   utils.track
