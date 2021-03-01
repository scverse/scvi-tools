=========
Developer
=========


Import scvi-tools as::

   import scvi


.. currentmodule:: scvi

Architectures
-------------

.. currentmodule:: scvi
   nn.FCLayers
   nn.Encoder
   nn.Decoder
   nn.one_hot
   nn.auto_move_data


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


.. autosummary::
    :toctree: reference/
    :nosignatures:

    model.base.BaseModelClass
    model.base.VAEMixin
    model.base.RNASeqMixin
    model.base.ArchesMixin


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


.. autosummary::
   :toctree: reference/
   :nosignatures:

   modules.base.LossRecorder
   modules.base.BaseModuleClass
   modules.base.PyroBaseModuleClass

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
