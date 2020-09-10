====
Core
====

.. warning::

    The ``scvi.core`` top-level module is PRIVATE. We document here for contributors. Please use the top-level module ``scvi.models``.

.. currentmodule:: scvi

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

Posteriors
~~~~~~~~~~

.. autosummary::
   :toctree: reference/

   core.posteriors.Posterior
   core.posteriors.TotalPosterior
   core.posteriors.AnnotationPosterior

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

Miscellaneous
~~~~~~~~~~~~~

.. autosummary::
   :toctree: reference/

   core.models.differential.DifferentialComputation
