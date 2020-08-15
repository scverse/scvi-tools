=======
Classes
=======

.. warning::

    The ``scvi.models._modules`` and ``scvi.inference`` modules are PRIVATE. We document here for completeness. Please use the top-level module ``scvi.models``.

.. currentmodule:: scvi

Modules
~~~~~~~

.. autosummary::
   :toctree: reference/

   models._modules.vae.VAE
   models._modules.vae.LDVAE
   models._modules.totalvae.TOTALVAE

Posteriors
~~~~~~~~~~

.. autosummary::
   :toctree: reference/

   inference.posterior.Posterior
   inference.total_inference.TotalPosterior

Trainers
~~~~~~~~

.. autosummary::
   :toctree: reference/

   inference.trainer.Trainer
   inference.inference.UnsupervisedTrainer
   inference.total_inference.TotalTrainer
   inference.trainer.EarlyStopping
