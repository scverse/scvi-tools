import copy
import logging
from typing import Union

import anndata
from numpy import ceil

from .trainer import Trainer

logger = logging.getLogger(__name__)


class UnsupervisedTrainer(Trainer):
    """
    Class for unsupervised training of an autoencoder.

    Parameters
    ----------
    model
        A model instance from class ``VAE``, ``VAEC``, ``SCANVI``, ``AutoZIVAE``
    adata
        A registered AnnData object
    train_size
        The train size, a float between 0 and 1 representing proportion of dataset to use for training
        to use Default: ``0.9``.
    test_size
        The test size,  a float between 0 and 1 representing proportion of dataset to use for testing
        to use Default: ``None``, which is equivalent to data not in the train set. If ``train_size`` and ``test_size``
        do not add to 1 then the remaining samples are added to a ``validation_set``.
    **kwargs
        Other keywords arguments from the general Trainer class.

    Other Parameters
    ----------------
    n_epochs_kl_warmup
        Number of epochs for linear warmup of KL(q(z|x)||p(z)) term. After `n_epochs_kl_warmup`,
        the training objective is the ELBO. This might be used to prevent inactivity of latent units, and/or to
        improve clustering of latent space, as a long warmup turns the model into something more of an autoencoder.
        Be aware that large datasets should avoid this mode and rely on n_iter_kl_warmup. If this parameter is not
        None, then it overrides any choice of `n_iter_kl_warmup`.
    n_iter_kl_warmup
        Number of iterations for warmup (useful for bigger datasets)
        int(128*5000/400) is a good default value.
    normalize_loss
        A boolean determining whether the loss is divided by the total number of samples used for
        training. In particular, when the global KL divergence is equal to 0 and the division is performed, the loss
        for a minibatchis is equal to the average of reconstruction losses and KL divergences on the minibatch.
        Default: ``None``, which is equivalent to setting False when the model is an instance from class
        ``AutoZIVAE`` and True otherwise.

    Examples
    --------
    >>> gene_dataset = CortexDataset()
    >>> vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
    ... n_labels=gene_dataset.n_labels)

    >>> infer = VariationalInference(gene_dataset, vae, train_size=0.5)
    >>> infer.train(n_epochs=20, lr=1e-3)

    Notes
    -----
    Two parameters can help control the training KL annealing
    If your applications rely on the posterior quality,
    (i.e. differential expression, batch effect removal), ensure the number of total
    epochs (or iterations) exceed the number of epochs (or iterations) used for KL warmup
    """

    default_metrics_to_monitor = ["elbo"]

    def __init__(
        self,
        model,
        adata: anndata.AnnData,
        train_size: Union[int, float] = 0.9,
        test_size: Union[int, float] = None,
        n_iter_kl_warmup: Union[int, None] = None,
        n_epochs_kl_warmup: Union[int, None] = 400,
        normalize_loss: bool = None,
        **kwargs
    ):
        train_size = float(train_size)
        if train_size > 1.0 or train_size <= 0.0:
            raise ValueError(
                "train_size needs to be greater than 0 and less than or equal to 1"
            )
        super().__init__(model, adata, **kwargs)

        # Set up number of warmup iterations
        self.n_iter_kl_warmup = n_iter_kl_warmup
        self.n_epochs_kl_warmup = n_epochs_kl_warmup
        self.normalize_loss = (
            not (
                hasattr(self.model, "gene_likelihood")
                and self.model.gene_likelihood == "autozinb"
            )
            if normalize_loss is None
            else normalize_loss
        )

        # Total size of the dataset used for training
        # (e.g. training set in this class but testing set in AdapterTrainer).
        # It used to rescale minibatch losses (cf. eq. (8) in Kingma et al., Auto-Encoding Variational Bayes, ICLR 2014)
        self.n_samples = 1.0

        if type(self) is UnsupervisedTrainer:
            (
                self.train_set,
                self.test_set,
                self.validation_set,
            ) = self.train_test_validation(model, adata, train_size, test_size)
            self.train_set.to_monitor = ["elbo"]
            self.test_set.to_monitor = ["elbo"]
            self.validation_set.to_monitor = ["elbo"]
            self.n_samples = len(self.train_set.indices)

    @property
    def scvi_data_loaders_loop(self):
        return ["train_set"]

    def loss(self, tensors: dict, feed_labels: bool = True):
        # The next lines should not be modified, because scanVI's trainer inherits
        # from this class and should NOT include label information to compute the ELBO by default
        loss_kwargs = dict(kl_weight=self.kl_weight, normalize_loss=self.normalize_loss)
        _, losses = self.model(tensors, loss_kwargs=loss_kwargs)
        loss = losses["loss"]
        return loss

    @property
    def kl_weight(self):
        epoch_criterion = self.n_epochs_kl_warmup is not None
        iter_criterion = self.n_iter_kl_warmup is not None
        if epoch_criterion:
            kl_weight = min(1.0, self.epoch / self.n_epochs_kl_warmup)
        elif iter_criterion:
            kl_weight = min(1.0, self.n_iter / self.n_iter_kl_warmup)
        else:
            kl_weight = 1.0
        return kl_weight

    def on_training_begin(self):
        epoch_criterion = self.n_epochs_kl_warmup is not None
        iter_criterion = self.n_iter_kl_warmup is not None
        if epoch_criterion:
            log_message = "KL warmup for {} epochs".format(self.n_epochs_kl_warmup)
            if self.n_epochs_kl_warmup > self.n_epochs:
                logger.info(
                    "KL warmup phase exceeds overall training phase"
                    "If your applications rely on the posterior quality, "
                    "consider training for more epochs or reducing the kl warmup."
                )
        elif iter_criterion:
            log_message = "KL warmup for {} iterations".format(self.n_iter_kl_warmup)
            n_iter_per_epochs_approx = ceil(self.adata.shape[0] / self.batch_size)
            n_total_iter_approx = self.n_epochs * n_iter_per_epochs_approx
            if self.n_iter_kl_warmup > n_total_iter_approx:
                logger.info(
                    "KL warmup phase may exceed overall training phase."
                    "If your applications rely on posterior quality, "
                    "consider training for more epochs or reducing the kl warmup."
                )
        else:
            log_message = "Training without KL warmup"
        logger.info(log_message)

    def on_training_end(self):
        if self.kl_weight < 0.99:
            logger.info(
                "Training is still in warming up phase. "
                "If your applications rely on the posterior quality, "
                "consider training for more epochs or reducing the kl warmup."
            )
        logger.info(
            "Training time:  %i s. / %i epochs"
            % (int(self.training_time), self.n_epochs)
        )


class AdapterTrainer(UnsupervisedTrainer):
    def __init__(self, model, adata, test_data_loader, frequency=5):
        super().__init__(model, adata, frequency=frequency)
        self.test_set = test_data_loader
        self.test_set.to_monitor = ["elbo"]
        self.params = list(self.model.z_encoder.parameters()) + list(
            self.model.l_encoder.parameters()
        )
        self.z_encoder_state = copy.deepcopy(model.z_encoder.state_dict())
        self.l_encoder_state = copy.deepcopy(model.l_encoder.state_dict())
        self.n_scale = len(self.test_set.indices)

    @property
    def scvi_data_loaders_loop(self):
        return ["test_set"]

    def train(self, n_path=10, n_epochs=50, **kwargs):
        for i in range(n_path):
            # Re-initialize to create new path
            self.model.z_encoder.load_state_dict(self.z_encoder_state)
            self.model.l_encoder.load_state_dict(self.l_encoder_state)
            super().train(n_epochs, params=self.params, **kwargs)

        return min(self.history["elbo_test_set"])
