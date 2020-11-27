====
Core
====

.. warning::

    The ``scvi.core`` top-level module is PRIVATE. We document here for contributors. Please use the top-level module ``scvi.model``.

.. currentmodule:: scvi

Distributions
~~~~~~~~~~~~

.. autosummary::
   :toctree: reference/

   core.distributions.NegativeBinomial
   core.distributions.NegativeBinomialMixture
   core.distributions.ZeroInflatedNegativeBinomial

Modules
~~~~~~~

.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst

   core.modules.VAE
   core.modules.LDVAE
   core.modules.TOTALVAE
   core.modules.SCANVAE
   core.modules.JVAE
   core.modules.AutoZIVAE
   core.modules.Classifier

Data Loaders
~~~~~~~~~~~~

.. autosummary::
   :toctree: reference/

   core.data_loaders.ScviDataLoader

Lightning
~~~~~~~~~

PyTorch lightning is used to train our modules. Tasks define train/test/val optimization
steps for modules like `TOTALVAE`, `SCANVAE`, etc.

.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst

   core.lightning.VAETask
   core.lightning.SemiSupervisedTask
   core.lightning.AdversarialTask
   core.lightning.Trainer

Utilities
~~~~~~~~~

.. autosummary::
   :toctree: reference/

   core.utils.DifferentialComputation
