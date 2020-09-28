====
Core
====

.. warning::

    The ``scvi.core`` top-level module is PRIVATE. We document here for contributors. Please use the top-level module ``scvi.models``.

.. currentmodule:: scvi

Distributions
~~~~~~~~~~~~

.. autosummary::
   :toctree: reference/

   core.distributions.NegativeBinomial
   core.distributions.NegativeBinomialMixture
   core.distributions.ZeroInflatedNegativeBinomial

Models
~~~~~~

.. autosummary::
    :toctree: reference/

    core.models.BaseModelClass
    core.models.VAEMixin
    core.models.RNASeqMixin

Modules
~~~~~~~

.. autosummary::
   :toctree: reference/

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
   core.data_loaders.TotalDataLoader
   core.data_loaders.AnnotationDataLoader

Trainers
~~~~~~~~

.. autosummary::
   :toctree: reference/

   core.trainers.UnsupervisedTrainer
   core.trainers.TotalTrainer
   core.trainers.SemiSupervisedTrainer
   core.trainers.ClassifierTrainer
   core.trainers.trainer.Trainer
   core.trainers.trainer.EarlyStopping

Utilities
~~~~~~~~~~~~~

.. autosummary::
   :toctree: reference/

   core.utils.DifferentialComputation
