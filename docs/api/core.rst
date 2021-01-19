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
   compose.SCVILoss
   compose.AbstractVAE

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

   dataloaders.ScviDataLoader
   dataloaders.ConcatDataLoader

Lightning
~~~~~~~~~

PyTorch lightning is used to train our modules. Tasks define train/test/val optimization
steps for modules like `TOTALVAE`, `SCANVAE`, etc.

.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst

   lightning.VAETask
   lightning.SemiSupervisedTask
   lightning.AdversarialTask
   lightning.Trainer

Utilities
~~~~~~~~~

.. autosummary::
   :toctree: reference/

   utils.DifferentialComputation
