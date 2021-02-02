====
Core
====

.. currentmodule:: scvi

Compose
~~~~~~~
.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst

   compose.FCLayers
   compose.Encoder
   compose.LossRecorder
   compose.BaseModuleClass

.. autosummary::
   :toctree: reference/

   compose.one_hot
   compose.auto_move_data

Distributions
~~~~~~~~~~~~

.. autosummary::
   :toctree: reference/

   distributions.NegativeBinomial
   distributions.NegativeBinomialMixture
   distributions.ZeroInflatedNegativeBinomial


Model (Base)
~~~~~~~~~~~~

.. autosummary::
    :toctree: reference/

    model.base.BaseModelClass
    model.base.VAEMixin
    model.base.RNASeqMixin
    model.base.ArchesMixin

Modules
~~~~~~~

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

.. autosummary::
   :toctree: reference/

   dataloaders.AnnDataLoader
   dataloaders.ConcatDataLoader

Lightning
~~~~~~~~~

PyTorch lightning is used to train our modules. TrainingPlans define train/test/val optimization
steps for modules like `TOTALVAE`, `SCANVAE`, etc.

.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst

   lightning.TrainingPlan
   lightning.SemiSupervisedTrainingPlan
   lightning.AdversarialTrainingPlan
   lightning.Trainer

Utilities
~~~~~~~~~

.. autosummary::
   :toctree: reference/

   utils.DifferentialComputation
