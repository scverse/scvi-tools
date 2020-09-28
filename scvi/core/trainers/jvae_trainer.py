import logging
from functools import partial
from itertools import cycle
from typing import Iterable, List, Tuple, Union

import anndata
import numpy as np
import torch
from torch import nn

from scvi import _CONSTANTS
from scvi.core._log_likelihood import compute_elbo
from scvi.core.data_loaders import ScviDataLoader
from scvi.core.trainers.trainer import Trainer

logger = logging.getLogger(__name__)


class JvaeDataLoader(ScviDataLoader):
    def __init__(self, *args, mode=0, **kwargs):
        super().__init__(*args, **kwargs)
        self.mode = mode

    def elbo(self) -> float:
        elbo = compute_elbo(self.model, self, mode=self.mode)
        logger.debug("ELBO : %.4f" % elbo)
        return elbo


class JVAETrainer(Trainer):
    """The trainer class for the unsupervised training of JVAE.

    Parameters
    ----------
    model
        A model instance from class ``JVAE``
    discriminator
        A model instance of a classifier (with logit output)
    dataset_list
        list of ScviDataset instance like ``[CortexDataset(), SmfishDataset()]``
    train_size
        Train-test split ratio in (0,1) to split cells
    kappa
        float to weight the discriminator loss
    n_epochs_kl_warmup
        Number of epochs for linear warmup of KL(q(z|x)||p(z)) term. After `n_epochs_kl_warmup`,
        the training objective is the ELBO. This might be used to prevent inactivity of latent units, and/or to
        improve clustering of latent space, as a long warmup turns the model into something more of an autoencoder.
    **kwargs
        Other keywords arguments from the general Trainer class.
    """

    default_metrics_to_monitor = ["elbo"]

    def __init__(
        self,
        model: nn.Module,
        discriminator: nn.Module,
        dataset_list: List[anndata.AnnData],
        train_size: float = 0.9,
        use_cuda: bool = True,
        kappa: float = 1.0,
        n_epochs_kl_warmup: int = 400,
        **kwargs
    ):
        train_size = float(train_size)
        if train_size > 1.0 or train_size <= 0.0:
            raise ValueError(
                "train_size needs to be greater than 0 and less than or equal to 1"
            )

        super().__init__(model, dataset_list[0], use_cuda=use_cuda, **kwargs)
        self.n_epochs_kl_warmup = n_epochs_kl_warmup
        self.kappa = kappa
        self.all_dataset = [
            self.create_scvi_dl(adata=gd, type_class=partial(JvaeDataLoader, mode=i))
            for i, gd in enumerate(dataset_list)
        ]
        self.n_dataset = len(self.all_dataset)
        self.all_train, self.all_test, self.all_validation = list(
            zip(
                *[
                    self.train_test_validation(
                        model,
                        gd,
                        train_size,
                        type_class=partial(JvaeDataLoader, mode=i),
                    )
                    for i, gd in enumerate(dataset_list)
                ]
            )
        )
        for i, d in enumerate(self.all_train):
            self.register_data_loader("train_%d" % i, d)
            d.to_monitor = ["elbo"]

        for i, d in enumerate(self.all_test):
            self.register_data_loader("test_%d" % i, d)
            d.to_monitor = ["elbo"]

        self.discriminator = discriminator
        if self.use_cuda:
            self.discriminator.cuda()

        self.kl_weight = None
        self.compute_metrics_time = None
        self.n_epochs = None

        self.track_disc = []

    def on_epoch_begin(self):
        if self.n_epochs_kl_warmup is not None:
            self.kl_weight = min(1, self.epoch / self.n_epochs_kl_warmup)
        else:
            self.kl_weight = 1.0

    @property
    def scvi_data_loaders_loop(self):
        return ["train_%d" % i for i in range(self.n_dataset)]

    def data_loaders_loop(self):
        scdls = [self._scvi_data_loaders[name] for name in self.scvi_data_loaders_loop]
        # find the largest dataset to cycle over the others
        largest = np.argmax([scdl.dataset.n_cells for scdl in scdls])

        data_loaders = [
            scdl if i == largest else cycle(scdl) for i, scdl in enumerate(scdls)
        ]

        return zip(*data_loaders)

    def on_training_loop(self, tensors_dict):

        if self.train_discriminator:
            latent_tensors = []
            for (i, tensors) in enumerate(tensors_dict):
                data = tensors[_CONSTANTS.X_KEY]
                z = self.model.sample_from_posterior_z(
                    data, mode=i, deterministic=False
                )
                latent_tensors.append(z)

            # Train discriminator
            d_loss = self.loss_discriminator([t.detach() for t in latent_tensors], True)
            d_loss *= self.kappa
            self.d_optimizer.zero_grad()
            d_loss.backward()
            self.d_optimizer.step()

            # Train generative model to fool discriminator
            fool_loss = self.loss_discriminator(latent_tensors, False)
            fool_loss *= self.kappa
            self.optimizer.zero_grad()
            fool_loss.backward()
            self.optimizer.step()

        # Train generative model
        self.current_loss = g_loss = self.loss(tensors_dict)
        self.optimizer.zero_grad()
        g_loss.backward()
        self.optimizer.step()

    def training_extras_init(self, lr_d=1e-3, eps=0.01):
        self.discriminator.train()

        d_params = filter(lambda p: p.requires_grad, self.discriminator.parameters())
        self.d_optimizer = torch.optim.Adam(d_params, lr=lr_d, eps=eps)
        self.train_discriminator = self.n_dataset > 1 and self.kappa > 0

    def training_extras_end(self):
        self.discriminator.eval()

    def loss_discriminator(
        self,
        latent_tensors: List[torch.Tensor],
        predict_true_class: bool = True,
        return_details: bool = False,
    ) -> Union[List[torch.Tensor], torch.Tensor]:
        """Compute the loss of the discriminator (either for the true labels or the fool labels)

        Parameters
        ----------
        latent_tensors
            Tensors for each dataset of the latent space
        predict_true_class
            Specify if the loss aims at minimizing the accuracy or the mixing
        return_details
            Boolean used to inspect the loss values, return detailed loss for each dataset
        latent_tensors: List[torch.Tensor]

        Returns
        -------
        type
            scalar loss if return_details is False, else list of scalar losses for each dataset

        """
        n_classes = self.n_dataset
        losses = []
        for i, z in enumerate(latent_tensors):
            cls_logits = nn.LogSoftmax(dim=1)(self.discriminator(z))

            if predict_true_class:
                cls_target = torch.zeros(
                    n_classes, dtype=torch.float32, device=z.device
                )
                cls_target[i] = 1.0
            else:
                cls_target = torch.ones(
                    n_classes, dtype=torch.float32, device=z.device
                ) / (n_classes - 1)
                cls_target[i] = 0.0

            l_soft = cls_logits * cls_target
            cls_loss = -l_soft.sum(dim=1).mean()
            losses.append(cls_loss)

        if return_details:
            return losses

        total_loss = torch.stack(losses).sum()
        return total_loss

    def loss(
        self, tensors: Iterable[torch.Tensor], return_details: bool = False
    ) -> Union[torch.Tensor, Tuple[List[torch.Tensor], List[torch.Tensor]]]:
        """Compute the loss of vae (reconstruction + kl_divergence)

        Parameters
        ----------
        tensors
            Tensors of observations for each dataset
        return_details
            Boolean used to inspect the loss values, return detailed loss for each dataset

        Returns
        -------
        type
            scalar loss if return_details is False, else tuple (reconstruction_loss, kl_loss)


        """
        reconstruction_losses = []
        kl_divergences = []
        losses = []
        total_batch_size = 0
        for i, data in enumerate(tensors):
            sample_batch, l_mean, l_var, batch_index, labels, *_ = self._unpack_tensors(
                data
            )
            reconstruction_loss, kl_divergence, _ = self.model(
                sample_batch, l_mean, l_var, batch_index, mode=i
            )
            loss = torch.mean(
                reconstruction_loss + self.kl_weight * kl_divergence
            ) * sample_batch.size(0)
            total_batch_size += sample_batch.size(0)
            losses.append(loss)
            if return_details:
                reconstruction_losses.append(reconstruction_loss.mean())
                kl_divergences.append(kl_divergence.mean())

        if return_details:
            return reconstruction_losses, kl_divergences

        averaged_loss = torch.stack(losses).sum() / total_batch_size
        return averaged_loss

    def _unpack_tensors(self, tensors):
        x = tensors[_CONSTANTS.X_KEY].squeeze_(0)
        local_l_mean = tensors[_CONSTANTS.LOCAL_L_MEAN_KEY].squeeze_(0)
        local_l_var = tensors[_CONSTANTS.LOCAL_L_VAR_KEY].squeeze_(0)
        batch_index = tensors[_CONSTANTS.BATCH_KEY].squeeze_(0)
        y = tensors[_CONSTANTS.LABELS_KEY].squeeze_(0)
        return x, local_l_mean, local_l_var, batch_index, y

    def get_discriminator_confusion(self) -> np.ndarray:
        """A good mixing should lead to a uniform matrix."""
        confusion = []
        for i, scdl in enumerate(self.all_dataset):

            indices = np.arange(scdl.dataset.n_cells)
            data = scdl.dataset[indices][_CONSTANTS.X_KEY]
            data = torch.from_numpy(data)
            if self.use_cuda:
                data = data.cuda()

            z = self.model.sample_from_posterior_z(data, mode=i, deterministic=True)
            cls_z = nn.Softmax(dim=1)(self.discriminator(z)).detach()

            cls_z = cls_z.cpu().numpy()

            row = cls_z.mean(axis=0)
            confusion.append(row)
        return np.array(confusion)

    def get_loss_magnitude(
        self, one_sample: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the different losses of the model separately. Useful to inspect and compare their magnitude.

        Parameters
        ----------
        one_sample
            Use only one batch to estimate the loss, can be much faster/less exact on big datasets
        """
        total_reconstruction = np.zeros(self.n_dataset)
        total_kl_divergence = np.zeros(self.n_dataset)
        total_discriminator = np.zeros(self.n_dataset)

        for tensors_dict in self.data_loaders_loop():
            reconstruction_losses, kl_divergences = self.loss(
                tensors_dict, return_details=True
            )
            discriminator_losses = self.loss_discriminator(
                [
                    self.model.sample_from_posterior_z(
                        data[_CONSTANTS.X_KEY], mode=i, deterministic=False
                    )
                    for (i, data) in enumerate(tensors_dict)
                ],
                return_details=True,
            )

            for i in range(self.n_dataset):
                total_reconstruction[i] += reconstruction_losses[i].item()
                total_kl_divergence[i] += kl_divergences[i].item()
                total_discriminator[i] += discriminator_losses[i].item()
            if one_sample:
                break

        return total_reconstruction, total_kl_divergence, total_discriminator
