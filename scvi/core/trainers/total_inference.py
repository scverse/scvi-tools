import logging
from typing import Union

import anndata
import torch

from scvi import _CONSTANTS
from scvi.core.data_loaders import TotalDataLoader
from scvi.core.modules import TOTALVAE, Classifier
from scvi.core.modules.utils import one_hot
from scvi.core.trainers import UnsupervisedTrainer

logger = logging.getLogger(__name__)


default_early_stopping_kwargs = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "patience": 45,
    "threshold": 0,
    "reduce_lr_on_plateau": True,
    "lr_patience": 30,
    "lr_factor": 0.6,
    "scvi_data_loader_class": TotalDataLoader,
}


def _unpack_tensors(tensors):
    x = tensors[_CONSTANTS.X_KEY]
    local_l_mean = tensors[_CONSTANTS.LOCAL_L_MEAN_KEY]
    local_l_var = tensors[_CONSTANTS.LOCAL_L_VAR_KEY]
    batch_index = tensors[_CONSTANTS.BATCH_KEY]
    labels = tensors[_CONSTANTS.LABELS_KEY]
    y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
    return x, local_l_mean, local_l_var, batch_index, labels, y


class TotalTrainer(UnsupervisedTrainer):
    """
    Unsupervised training for totalVI using variational inference.

    Parameters
    ----------
    model
        A model instance from class ``TOTALVAE``
    adata
        A registered AnnData object
    train_size
        The train size, a float between 0 and 1 representing proportion of dataset to use for training
        to use Default: ``0.90``.
    test_size
        The test size, a float between 0 and 1 representing proportion of dataset to use for testing
        to use Default: ``0.10``. Note that if train and test do not add to 1 the remainder is placed in a validation set
    pro_recons_weight
        Scaling factor on the reconstruction loss for proteins. Default: ``1.0``.
    n_epochs_kl_warmup
        Number of epochs for annealing the KL terms for `z` and `mu` of the ELBO (from 0 to 1). If None, no warmup performed, unless
        `n_iter_kl_warmup` is set.
    n_iter_kl_warmup
        Number of minibatches for annealing the KL terms for `z` and `mu` of the ELBO (from 0 to 1). If set to "auto", the number
        of iterations is equal to 75% of the number of cells. `n_epochs_kl_warmup` takes precedence if it is not None. If both are None, then
        no warmup is performed.
    discriminator
        Classifier used for adversarial training scheme
    use_adversarial_loss
        Whether to use adversarial classifier to improve mixing
    kappa
        Scaling factor for adversarial loss. If None, follow inverse of kl warmup schedule.
    early_stopping_kwargs
        Keyword args for early stopping. If "auto", use totalVI defaults. If None, disable early stopping.
    """

    default_metrics_to_monitor = ["elbo"]

    def __init__(
        self,
        model: TOTALVAE,
        dataset: anndata.AnnData,
        train_size: float = 0.90,
        test_size: float = 0.10,
        pro_recons_weight: float = 1.0,
        n_epochs_kl_warmup: int = None,
        n_iter_kl_warmup: Union[str, int] = "auto",
        discriminator: Classifier = None,
        use_adversarial_loss: bool = False,
        kappa: float = None,
        early_stopping_kwargs: Union[dict, str, None] = "auto",
        **kwargs,
    ):
        train_size = float(train_size)
        if train_size > 1.0 or train_size <= 0.0:
            raise ValueError(
                "train_size needs to be greater than 0 and less than or equal to 1"
            )
        self.n_genes = model.n_input_genes
        self.n_proteins = model.n_input_proteins
        self.use_adversarial_loss = use_adversarial_loss
        self.kappa = kappa
        self.pro_recons_weight = pro_recons_weight

        if early_stopping_kwargs == "auto":
            early_stopping_kwargs = default_early_stopping_kwargs

        super().__init__(
            model,
            dataset,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            n_iter_kl_warmup=0.75 * len(dataset)
            if n_iter_kl_warmup == "auto"
            else n_iter_kl_warmup,
            early_stopping_kwargs=early_stopping_kwargs,
            **kwargs,
        )

        if use_adversarial_loss is True and discriminator is None:
            discriminator = Classifier(
                n_input=self.model.n_latent,
                n_hidden=32,
                n_labels=self.adata.uns["_scvi"]["summary_stats"]["n_batch"],
                n_layers=2,
                logits=True,
            )

        self.discriminator = discriminator
        if self.use_cuda and self.discriminator is not None:
            self.discriminator.cuda()

        if isinstance(self, TotalTrainer):
            (
                self.train_set,
                self.test_set,
                self.validation_set,
            ) = self.train_test_validation(
                model, dataset, train_size, test_size, type_class=TotalDataLoader
            )
            self.train_set.to_monitor = []
            self.test_set.to_monitor = ["elbo"]
            self.validation_set.to_monitor = ["elbo"]

    def loss(self, tensors):
        (
            sample_batch_x,
            local_l_mean,
            local_l_var,
            batch_index,
            label,
            sample_batch_y,
        ) = _unpack_tensors(tensors)
        (
            reconst_loss_gene,
            reconst_loss_protein,
            kl_div_z,
            kl_div_l_gene,
            kl_div_back_pro,
        ) = self.model(
            sample_batch_x,
            sample_batch_y,
            local_l_mean,
            local_l_var,
            batch_index,
            label,
        )

        loss = torch.mean(
            reconst_loss_gene
            + self.pro_recons_weight * reconst_loss_protein
            + self.kl_weight * kl_div_z
            + kl_div_l_gene
            + self.kl_weight * kl_div_back_pro
        )

        return loss

    def loss_discriminator(
        self, z, batch_index, predict_true_class=True, return_details=True
    ):

        n_classes = self.adata.uns["_scvi"]["summary_stats"]["n_batch"]
        cls_logits = torch.nn.LogSoftmax(dim=1)(self.discriminator(z))

        if predict_true_class:
            cls_target = one_hot(batch_index, n_classes)
        else:
            one_hot_batch = one_hot(batch_index, n_classes)
            cls_target = torch.zeros_like(one_hot_batch)
            # place zeroes where true label is
            cls_target.masked_scatter_(
                ~one_hot_batch.bool(), torch.ones_like(one_hot_batch) / (n_classes - 1)
            )

        l_soft = cls_logits * cls_target
        loss = -l_soft.sum(dim=1).mean()

        return loss

    def _get_z(self, tensors):
        (
            sample_batch_x,
            local_l_mean,
            local_l_var,
            batch_index,
            label,
            sample_batch_y,
        ) = _unpack_tensors(tensors)

        z = self.model.sample_from_posterior_z(
            sample_batch_x, sample_batch_y, batch_index, give_mean=False
        )

        return z

    def train(self, n_epochs=500, lr=4e-3, eps=0.01, params=None, max_grad_value=None):
        self.max_grad_value = max_grad_value

        super().train(n_epochs=n_epochs, lr=lr, eps=eps, params=params)

    def on_training_loop(self, tensors_dict):
        if self.use_adversarial_loss:
            if self.kappa is None:
                kappa = 1 - self.kl_weight
            else:
                kappa = self.kappa
            batch_index = tensors_dict[0][_CONSTANTS.BATCH_KEY]
            if kappa > 0:
                z = self._get_z(*tensors_dict)
                # Train discriminator
                d_loss = self.loss_discriminator(z.detach(), batch_index, True)
                d_loss *= kappa
                self.d_optimizer.zero_grad()
                d_loss.backward()
                self.d_optimizer.step()

                # Train generative model to fool discriminator
                fool_loss = self.loss_discriminator(z, batch_index, False)
                fool_loss *= kappa

            # Train generative model
            self.optimizer.zero_grad()
            self.current_loss = loss = self.loss(*tensors_dict)
            if kappa > 0:
                (loss + fool_loss).backward()
            else:
                loss.backward()
            if self.max_grad_value is not None:
                torch.nn.utils.clip_grad_norm_(
                    self.optimizer.param_groups[0]["params"],
                    self.max_grad_value,
                    norm_type="inf",
                )
            self.optimizer.step()

        else:
            self.current_loss = loss = self.loss(*tensors_dict)
            self.optimizer.zero_grad()
            loss.backward()
            if self.max_grad_value is not None:
                torch.nn.utils.clip_grad_norm_(
                    self.optimizer.param_groups[0]["params"],
                    self.max_grad_value,
                    norm_type="inf",
                )
            self.optimizer.step()

    def training_extras_init(self, lr_d=1e-3, eps=0.01):
        if self.discriminator is not None:
            self.discriminator.train()

            d_params = filter(
                lambda p: p.requires_grad, self.discriminator.parameters()
            )
            self.d_optimizer = torch.optim.Adam(d_params, lr=lr_d, eps=eps)

    def training_extras_end(self):
        if self.discriminator is not None:
            self.discriminator.eval()
