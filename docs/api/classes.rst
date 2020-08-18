=======
Classes
=======

.. warning::

    The ``scvi.core`` top-level module is PRIVATE. We document here for developers. Please use the top-level module ``scvi.models``.

.. currentmodule:: scvi

Modules
~~~~~~~

.. autosummary::
   :toctree: reference/

   core.models.VAE
   core.models.LDVAE
   core.models.TOTALVAE

Posteriors
~~~~~~~~~~

.. autosummary::
   :toctree: reference/

   core.posteriors.Posterior
   core.posteriors.TotalPosterior

Trainers
~~~~~~~~

.. autosummary::
   :toctree: reference/

   core.trainers.UnsupervisedTrainer
   core.trainers.TotalTrainer
   core.trainers.trainer.Trainer
   core.trainers.trainer.EarlyStopping
