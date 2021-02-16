=========
Developer
=========


Import scvi-tools as::

   import scvi


.. currentmodule:: scvi

Compose
~~~~~~~

.. currentmodule:: scvi

.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst

   compose.FCLayers
   compose.Encoder
   compose.LossRecorder
   compose.BaseModuleClass
   compose.PyroBaseModuleClass

.. autosummary::
   :toctree: reference/

   compose.one_hot
   compose.auto_move_data

Distributions
~~~~~~~~~~~~~

.. currentmodule:: scvi


.. autosummary::
   :toctree: reference/

   distributions.NegativeBinomial
   distributions.NegativeBinomialMixture
   distributions.ZeroInflatedNegativeBinomial


Model (Base)
~~~~~~~~~~~~

.. currentmodule:: scvi


.. autosummary::
    :toctree: reference/

    model.base.BaseModelClass
    model.base.VAEMixin
    model.base.RNASeqMixin
    model.base.ArchesMixin

Modules
~~~~~~~

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
~~~~~~~~~~~~

.. currentmodule:: scvi


.. autosummary::
   :toctree: reference/

   dataloaders.AnnDataLoader
   dataloaders.AnnTorchDataset
   dataloaders.ConcatDataLoader
   dataloaders.SemiSupervisedDataLoader

Lightning
~~~~~~~~~

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
~~~~~~~~~

.. currentmodule:: scvi


.. autosummary::
   :toctree: reference/

   utils.DifferentialComputation
