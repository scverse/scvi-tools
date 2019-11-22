import logging
import copy
import warnings
from typing import Union

import matplotlib.pyplot as plt
import torch
from numpy import ceil

from scvi.dataset import GeneExpressionDataset
from scvi.inference import Trainer

plt.switch_backend("agg")
logger = logging.getLogger(__name__)


class UnsupervisedTrainer(Trainer):
    r"""The VariationalInference class for the unsupervised training of an autoencoder.

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SCANVI``, ``AutoZIVAE``
        :gene_dataset: A gene_dataset instance like ``CortexDataset()``
        :train_size: The train size, either a float between 0 and 1 or an integer for the number of training samples
         to use Default: ``0.8``.
        :test_size: The test size, either a float between 0 and 1 or an integer for the number of training samples
         to use Default: ``None``, which is equivalent to data not in the train set. If ``train_size`` and ``test_size``
         do not add to 1 or the length of the dataset then the remaining samples are added to a ``validation_set``.
        :n_epochs_kl_warmup: Number of epochs for linear warmup of KL(q(z|x)||p(z)) term. After `n_epochs_kl_warmup`,
            the training objective is the ELBO. This might be used to prevent inactivity of latent units, and/or to
            improve clustering of latent space, as a long warmup turns the model into something more of an autoencoder.
        :normalize_loss: A boolean determining whether the loss is divided by the total number of samples used for
            training. In particular, when the global KL divergence is equal to 0 and the division is performed, the loss
            for a minibatchis is equal to the average of reconstruction losses and KL divergences on the minibatch.
            Default: ``None``, which is equivalent to setting False when the model is an instance from class
            ``AutoZIVAE`` and True otherwise.
        :\*\*kwargs: Other keywords arguments from the general Trainer class.

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
        ... n_labels=gene_dataset.n_labels)

        >>> infer = VariationalInference(gene_dataset, vae, train_size=0.5)
        >>> infer.train(n_epochs=20, lr=1e-3)
    """
    default_metrics_to_monitor = ["elbo"]

    def __init__(
        self,
        model,
        gene_dataset: GeneExpressionDataset,
        train_size: Union[int, float] = 0.8,
        test_size: Union[int, float] = None,
        n_iter_kl_warmup: int = 50,
        n_epochs_kl_warmup: int = None,
        normalize_loss: bool = None,
        **kwargs
    ):
        super().__init__(model, gene_dataset, **kwargs)

        # Set up number of warmup iterations
        n_reference_cells = min(10000, gene_dataset.nb_cells)
        self.n_iter_kl_warmup = None
        self.n_epochs_kl_warmup = None

        if n_iter_kl_warmup is not None:
            self.n_iter_kl_warmup = (
                1.0 * n_iter_kl_warmup * ceil(n_reference_cells / self.batch_size)
            )
        elif n_epochs_kl_warmup is not None:
            self.n_epochs_kl_warmup = n_epochs_kl_warmup
        self.n_epochs_kl_warmup = n_epochs_kl_warmup
        self.normalize_loss = (
            not (
                hasattr(self.model, "reconstruction_loss")
                and self.model.reconstruction_loss == "autozinb"
            )
            if normalize_loss is None
            else normalize_loss
        )

        # Total size of the dataset used for training
        # (e.g. training set in this class but testing set in AdapterTrainer).
        # It used to rescale minibatch losses (cf. eq. (8) in Kingma et al., Auto-Encoding Variational Bayes, iCLR 2013)
        self.n_samples = 1.0

        if type(self) is UnsupervisedTrainer:
            (
                self.train_set,
                self.test_set,
                self.validation_set,
            ) = self.train_test_validation(model, gene_dataset, train_size, test_size)
            self.train_set.to_monitor = ["elbo"]
            self.test_set.to_monitor = ["elbo"]
            self.validation_set.to_monitor = ["elbo"]
            self.n_samples = len(self.train_set.indices)

    @property
    def posteriors_loop(self):
        return ["train_set"]

    def loss(self, tensors):
        sample_batch, local_l_mean, local_l_var, batch_index, y = tensors
        reconst_loss, kl_divergence_local, kl_divergence_global = self.model(
            sample_batch, local_l_mean, local_l_var, batch_index, y
        )
        loss = (
            self.n_samples
            * torch.mean(reconst_loss + self.kl_weight * kl_divergence_local)
            + kl_divergence_global
        )
        if self.normalize_loss:
            loss = loss / self.n_samples
        return loss

    @property
    def kl_weight(self):
        epoch_criterium = self.n_epochs_kl_warmup is not None
        iter_criterium = self.n_iter_kl_warmup is not None
        if epoch_criterium:
            kl_weight = min(1.0, self.epoch / self.n_epochs_kl_warmup)
        elif iter_criterium:
            kl_weight = min(1.0, self.n_iter / self.n_iter_kl_warmup)
        else:
            kl_weight = 1.0
        return kl_weight

    def on_training_begin(self):
        epoch_criterium = self.n_epochs_kl_warmup is not None
        iter_criterium = self.n_iter_kl_warmup is not None
        if iter_criterium:
            log_message = "KL warmup for {} iterations".format(self.n_iter_kl_warmup)
        elif epoch_criterium:
            log_message = "KL warmup for {} epochs".format(self.n_epochs_kl_warmup)
            if self.n_epochs_kl_warmup > self.n_epochs:
                raise ValueError("KL warmup phase exceeds overall training phase")
        else:
            log_message = "Training without KL warmup"
        logger.info(log_message)

    def on_training_end(self):
        if self.kl_weight < 0.99:
            warnings.warn(
                "Training is still is warming up phase. Consider increasing train's n_epochs "
                "or reducing n_iter_kl_warmup/n_epochs_kl_warmup"
            )


class AdapterTrainer(UnsupervisedTrainer):
    def __init__(self, model, gene_dataset, posterior_test, frequency=5):
        super().__init__(model, gene_dataset, frequency=frequency)
        self.test_set = posterior_test
        self.test_set.to_monitor = ["elbo"]
        self.params = list(self.model.z_encoder.parameters()) + list(
            self.model.l_encoder.parameters()
        )
        self.z_encoder_state = copy.deepcopy(model.z_encoder.state_dict())
        self.l_encoder_state = copy.deepcopy(model.l_encoder.state_dict())
        self.n_scale = len(self.test_set.indices)

    @property
    def posteriors_loop(self):
        return ["test_set"]

    def train(self, n_path=10, n_epochs=50, **kwargs):
        for i in range(n_path):
            # Re-initialize to create new path
            self.model.z_encoder.load_state_dict(self.z_encoder_state)
            self.model.l_encoder.load_state_dict(self.l_encoder_state)
            super().train(n_epochs, params=self.params, **kwargs)

        return min(self.history["elbo_test_set"])
