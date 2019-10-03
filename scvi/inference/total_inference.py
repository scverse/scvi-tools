from typing import Optional
import logging
import torch
from torch.distributions import Poisson, Gamma, Bernoulli, Normal
from torch.utils.data import DataLoader
import numpy as np

from scvi.inference import Posterior
from . import UnsupervisedTrainer

from scvi.dataset import GeneExpressionDataset
from scvi.models import TOTALVI

logger = logging.getLogger(__name__)


class TotalPosterior(Posterior):
    r"""The functional data unit for totalVI. A `TotalPosterior` instance is instantiated with a model and
    a gene_dataset, and as well as additional arguments that for Pytorch's `DataLoader`. A subset of indices
    can be specified, for purposes such as splitting the data into train/test or labelled/unlabelled
    (for semi-supervised learning). Each trainer instance of the `TotalTrainer` class can therefore have multiple
    `TotalPosterior` instances to train a model. A `TotalPosterior` instance also comes with many methods or
    utilities for its corresponding data.


    :param model: A model instance from class ``TOTALVI``
    :param gene_dataset: A gene_dataset instance like ``CbmcDataset()``
    :param shuffle: Specifies if a `RandomSampler` or a `SequentialSampler` should be used
    :param indices: Specifies how the data should be split with regards to train/test or labelled/unlabelled
    :param use_cuda: Default: ``True``
    :param data_loader_kwarg: Keyword arguments to passed into the `DataLoader`

    Examples:

    Let us instantiate a `trainer`, with a gene_dataset and a model

        >>> gene_dataset = CbmcDataset()
        >>> totalvi = TOTALVI(gene_dataset.nb_genes, len(gene_dataset.protein_names),
        ... n_batch=gene_dataset.n_batches * False, n_labels=gene_dataset.n_labels, use_cuda=True)
        >>> trainer = TotalTrainer(vae, gene_dataset)
        >>> trainer.train(n_epochs=400)
    """

    def __init__(
        self,
        model: TOTALVI,
        gene_dataset: GeneExpressionDataset,
        shuffle: bool = False,
        indices: Optional[np.ndarray] = None,
        use_cuda: bool = True,
        data_loader_kwargs=dict(),
    ):

        super().__init__(
            model,
            gene_dataset,
            shuffle=shuffle,
            indices=indices,
            use_cuda=use_cuda,
            data_loader_kwargs=data_loader_kwargs,
        )
        # Add protein tensor as another tensor to be loaded
        self.data_loader_kwargs.update(
            {
                "collate_fn": gene_dataset.collate_fn_builder(
                    {"protein_expression": np.float32}
                )
            }
        )
        self.data_loader = DataLoader(gene_dataset, **self.data_loader_kwargs)

    def corrupted(self):
        return self.update(
            {
                "collate_fn": self.gene_dataset.collate_fn_builder(
                    {"protein_expression": np.float32}, corrupted=True
                )
            }
        )

    def uncorrupted(self):
        return self.update(
            {
                "collate_fn": self.gene_dataset.collate_fn_builder(
                    {"protein_expression": np.float32}
                )
            }
        )

    @torch.no_grad()
    def elbo(self):
        elbo = self.compute_elbo(self.model)
        return elbo

    elbo.mode = "min"

    @torch.no_grad()
    def reconstruction_error(self):
        ll_gene, ll_protein = self.compute_reconstruction_error(self.model)
        return ll_gene + ll_protein

    reconstruction_error.mode = "min"

    @torch.no_grad()
    def reconstruction_error_protein(self):
        ll_gene, ll_protein = self.compute_reconstruction_error(self.model)
        return ll_protein

    reconstruction_error_protein.mode = "min"

    @torch.no_grad()
    def reconstruction_error_separated(self):
        ll_gene, ll_protein = self.compute_reconstruction_error(self.model)
        return ll_gene, ll_protein

    @torch.no_grad()
    def marginal_ll(self, n_mc_samples=1000):
        raise NotImplementedError

    @torch.no_grad()
    def get_background_mean(self):
        background_mean = []
        for tensors in self:
            x, _, _, batch_index, labels, y = tensors
            background_mean += [
                np.array(
                    (
                        self.model.get_background_mean(
                            x, y, batch_index=batch_index, label=labels, n_samples=1
                        )
                    ).cpu()
                )
            ]
        return np.concatenate(background_mean)

    def compute_elbo(self, vae: TOTALVI, **kwargs):
        """ Computes the ELBO.

        The ELBO is the reconstruction error + the KL divergences
        between the variational distributions and the priors.
        It differs from the marginal log likelihood.
        Specifically, it is a lower bound on the marginal log likelihood
        plus a term that is constant with respect to the variational distribution.
        It still gives good insights on the modeling of the data, and is fast to compute.
        """
        # Iterate once over the posterior and computes the total log_likelihood
        elbo = 0
        for i_batch, tensors in enumerate(self):
            x, local_l_mean, local_l_var, batch_index, labels, y = tensors
            reconst_loss_gene, reconst_loss_protein, kl_divergence, back_kl = vae(
                x,
                y,
                local_l_mean,
                local_l_var,
                batch_index=batch_index,
                label=labels,
                **kwargs
            )
            elbo += torch.sum(
                reconst_loss_gene + reconst_loss_protein + kl_divergence + back_kl
            ).item()
        n_samples = len(self.indices)
        return elbo / n_samples

    def compute_reconstruction_error(self, vae: TOTALVI, **kwargs):
        r""" Computes log p(x/z), which is the reconstruction error .
            Differs from the marginal log likelihood, but still gives good
            insights on the modeling of the data, and is fast to compute

            This is really a helper function to self.ll, self.ll_protein, etc.
        """
        # Iterate once over the posterior and computes the total log_likelihood
        log_lkl_gene = 0
        log_lkl_protein = 0
        for i_batch, tensors in enumerate(self):
            x, local_l_mean, local_l_var, batch_index, labels, y = tensors
            reconst_loss_gene, reconst_loss_protein, kl_divergence, back_kl = vae(
                x,
                y,
                local_l_mean,
                local_l_var,
                batch_index=batch_index,
                label=labels,
                **kwargs
            )
            log_lkl_gene += torch.sum(reconst_loss_gene).item()
            log_lkl_protein += torch.sum(reconst_loss_protein).item()

        n_samples = len(self.indices)
        return log_lkl_gene / n_samples, log_lkl_protein / n_samples

    def compute_marginal_log_likelihood(
        self, n_samples_mc: int = 100, batch_size: int = 96
    ):
        """ Computes a biased estimator for log p(x, y), which is the marginal log likelihood.

        Despite its bias, the estimator still converges to the real value
        of log p(x, y) when n_samples_mc (for Monte Carlo) goes to infinity
        (a fairly high value like 100 should be enough)
        Due to the Monte Carlo sampling, this method is not as computationally efficient
        as computing only the reconstruction loss
        """
        # Uses MC sampling to compute a tighter lower bound on log p(x)
        log_lkl = 0
        for i_batch, tensors in enumerate(self.update({"batch_size": batch_size})):
            x, local_l_mean, local_l_var, batch_index, labels, y = tensors
            to_sum = torch.zeros(x.size()[0], n_samples_mc)

            for i in range(n_samples_mc):

                # Distribution parameters and sampled variables
                outputs = self.model.inference(x, y, batch_index, labels)
                qz_m = outputs["qz_m"]
                qz_v = outputs["qz_v"]
                ql_m = outputs["ql_m"]
                ql_v = outputs["ql_v"]
                px_rate = outputs["px_rate"]
                px_scale = outputs["px_scale"]
                px_r = outputs["px_r"]
                px_dropout = outputs["px_dropout"]
                library_gene = outputs["library_gene"]
                # really need not softmax transformed random variable
                z = outputs["untran_z"]
                log_back_mean = outputs["log_back_mean"]

                # Reconstruction Loss
                reconst_loss_gene, reconst_loss_protein = self.model.get_reconstruction_loss(
                    x, y, px_rate, px_r, px_dropout, px_scale
                )

                # Log-probabilities
                p_l_gene = (
                    Normal(local_l_mean, local_l_var.sqrt())
                    .log_prob(library_gene)
                    .sum(dim=-1)
                )
                p_z = self.model.z_prior.log_prob(z)
                if self.model.latent_distribution != "ln":
                    p_z = p_z.sum(dim=-1)
                p_mu_back = self.model.back_mean_prior.log_prob(log_back_mean).sum(
                    dim=-1
                )
                p_xy_zl = -(reconst_loss_gene + reconst_loss_protein)
                q_z_x = Normal(qz_m, qz_v.sqrt()).log_prob(z).sum(dim=-1)
                q_l_x = Normal(ql_m, ql_v.sqrt()).log_prob(library_gene).sum(dim=-1)
                q_mu_back = (
                    Normal(outputs["py_alpha"], outputs["py_beta"])
                    .log_prob(log_back_mean)
                    .sum(dim=-1)
                )
                to_sum[:, i] = (
                    p_z + p_l_gene + p_mu_back + p_xy_zl - q_z_x - q_l_x - q_mu_back
                )

            batch_log_lkl = torch.logsumexp(to_sum, dim=-1) - np.log(n_samples_mc)
            log_lkl += torch.sum(batch_log_lkl).item()

        n_samples = len(self.indices)
        # The minus sign is there because we actually look at the negative log likelihood
        return -log_lkl / n_samples

    @torch.no_grad()
    def get_latent(self, sample: bool = False):
        """
        Output posterior z mean or sample, batch index, and label
        :param sample: z mean or z sample
        :return: six np.ndarrays, latent, batch_indices, labels, library for each of gene, protein, background
        """
        latent = []
        batch_indices = []
        labels = []
        library_gene = []
        for tensors in self:
            x, local_l_mean, local_l_var, batch_index, label, y = tensors
            give_mean = not sample
            latent += [
                self.model.sample_from_posterior_z(
                    x, y, batch_index, give_mean=give_mean
                ).cpu()
            ]
            batch_indices += [batch_index.cpu()]
            labels += [label.cpu()]
            library_gene += [
                self.model.sample_from_posterior_l(
                    x, y, batch_index, give_mean=give_mean
                ).cpu()
            ]
        return (
            np.array(torch.cat(latent)),
            np.array(torch.cat(batch_indices)),
            np.array(torch.cat(labels)).ravel(),
            np.array(torch.cat(library_gene)).ravel(),
        )

    @torch.no_grad()
    def differential_expression_stats(self, M_sampling: int = 100):
        raise NotImplementedError

    @torch.no_grad()
    def get_harmonized_scale(self, fixed_batch: torch.Tensor):
        px_scales = []
        fixed_batch = float(fixed_batch)
        for tensors in self:
            x, local_l_mean, local_l_var, batch_index, label, y = tensors
            px_scales += [self.model.scale_from_z(x, y, fixed_batch).cpu()]
        return np.concatenate(px_scales)

    @torch.no_grad()
    def generate(
        self,
        n_samples: int = 100,
        genes: Optional[np.ndarray] = None,
        batch_size: int = 64,
    ):  # with n_samples>1 return original list/ otherwose sequential
        """
        Return original_values as y and generated as x (for posterior density visualization)
        :param n_samples:
        :param genes:
        :return:
        """
        original_list = []
        posterior_list = []
        for tensors in self.update({"batch_size": batch_size}):
            x, _, _, batch_index, labels, y = tensors
            with torch.no_grad():
                outputs = self.model.inference(
                    x, y, batch_index=batch_index, label=labels, n_samples=n_samples
                )
            px_dispersion = outputs["px_r"]
            px_rate = outputs["px_rate"]
            px_dropout = outputs["px_dropout"]

            pi = 1 / (1 + torch.exp(-px_dropout["protein"]))
            mixing_sample = Bernoulli(pi).sample()
            protein_rate = (
                px_rate["protein"] * (1 - mixing_sample)
                + px_rate["background"] * mixing_sample
            )
            px_rate = torch.cat((px_rate["gene"], protein_rate), dim=-1)
            if len(px_dispersion["gene"].size()) == 2:
                px_dispersion_gene = px_dispersion["gene"]
            else:
                px_dispersion_gene = torch.ones_like(x) * px_dispersion["gene"]
            if len(px_dispersion["protein"].size()) == 2:
                px_dispersion_protein = px_dispersion["protein"]
            else:
                px_dispersion_protein = torch.ones_like(y) * px_dispersion["protein"]

            px_dispersion = torch.cat(
                (px_dispersion_gene, px_dispersion_protein), dim=-1
            )

            # This gamma is really l*w using scVI manuscript notation
            p = px_rate / (px_rate + px_dispersion)
            r = px_dispersion
            l_train = Gamma(r, (1 - p) / p).sample()
            X = Poisson(l_train).sample().cpu().numpy()
            # """
            # In numpy (shape, scale) => (concentration, rate), with scale = p /(1 - p)
            # rate = (1 - p) / p  # = 1/scale # used in pytorch
            # """
            original_list += [np.array(torch.cat((x, y), dim=-1).cpu())]
            posterior_list += [X]

            if genes is not None:
                posterior_list[-1] = posterior_list[-1][
                    :, :, self.gene_dataset._gene_idx(genes)
                ]
                original_list[-1] = original_list[-1][
                    :, self.gene_dataset._gene_idx(genes)
                ]

            posterior_list[-1] = np.transpose(posterior_list[-1], (1, 2, 0))

        return (
            np.concatenate(posterior_list, axis=0),
            np.concatenate(original_list, axis=0),
        )

    @torch.no_grad()
    def get_sample_dropout(self, n_samples: int = 1, give_mean: bool = True):
        px_dropouts = []
        for tensors in self:
            x, _, _, batch_index, labels, y = tensors
            px_dropouts += [
                self.model.get_sample_dropout(
                    x, y, batch_index=batch_index, label=labels, n_samples=n_samples
                ).cpu()
            ]
        if n_samples > 1:
            # concatenate along batch dimension -> result shape = (samples, cells, features)
            px_dropouts = torch.cat(px_dropouts, dim=1)
            # (cells, features, samples)
            px_dropouts = px_dropouts.permute(1, 2, 0)
        else:
            px_dropouts = torch.cat(px_dropouts, dim=0)

        if give_mean is True and n_samples > 1:
            px_dropouts = torch.mean(px_dropouts, dim=-1)

        px_dropouts = px_dropouts.cpu().numpy()

        return px_dropouts

    @torch.no_grad()
    def get_normalized_denoised_expresssion(
        self, n_samples: int = 1, give_mean: bool = True
    ):
        """Returns the tensor denoised normalized gene and protein expression

        The protein expression columns are concatenated to the end of the gene expression

        :param n_samples: number of samples from posterior distribution
        :param give_mean: bool, whether to return samples along first axis or average over samples
        :rtype: :py:class:`np.ndarray`
        """

        px_scale_list = []
        for tensors in self:
            x, _, _, batch_index, labels, y = tensors
            outputs = self.model.inference(
                x, y, batch_index=batch_index, label=labels, n_samples=n_samples
            )
            px_scale = outputs["px_scale"]

            adjusted_px_scale = torch.cat(
                (px_scale["gene"], px_scale["protein"]), dim=-1
            )

            px_scale_list.append(adjusted_px_scale.cpu())

        if n_samples > 1:
            # concatenate along batch dimension -> result shape = (samples, cells, features)
            px_scale_list = torch.cat(px_scale_list, dim=1)
            # (cells, features, samples)
            px_scale_list = px_scale_list.permute(1, 2, 0)
        else:
            px_scale_list = torch.cat(px_scale_list, dim=0)

        if give_mean is True and n_samples > 1:
            px_scale_list = torch.mean(px_scale_list, dim=-1)

        px_scale_list = px_scale_list.cpu().numpy()

        return px_scale_list

    @torch.no_grad()
    def imputation(self, n_samples: int = 1):
        imputed_list = []
        for tensors in self:
            x, _, _, batch_index, labels, y = tensors
            px_rate = self.model.get_sample_rate(
                x, y, batch_index=batch_index, label=labels, n_samples=n_samples
            )
            imputed_list += [np.array(px_rate.cpu())]
        imputed_list = np.concatenate(imputed_list)
        return imputed_list.squeeze()

    @torch.no_grad()
    def get_sample_protein_rates(self, n_samples: int = 1, give_mean: bool = False):
        px_rate_list_foreground = []
        px_rate_list_background = []
        for tensors in self:
            x, _, _, batch_index, labels, y = tensors
            outputs = self.model.inference(
                x, y, batch_index=batch_index, label=labels, n_samples=n_samples
            )
            px_rate = outputs["px_rate"]
            px_rate_list_foreground.append(px_rate["protein"].cpu())
            px_rate_list_background.append(px_rate["background"].cpu())
        if n_samples > 1:
            # concatenate along batch dimension -> result shape = (samples, cells, features)
            px_rate_list_foreground = torch.cat(px_rate_list_foreground, dim=1)
            px_rate_list_background = torch.cat(px_rate_list_background, dim=1)
            # (cells, features, samples)
            px_rate_list_foreground = px_rate_list_foreground.permute(1, 2, 0)
            px_rate_list_background = px_rate_list_background.permute(1, 2, 0)
        else:
            px_rate_list_foreground = torch.cat(px_rate_list_foreground, dim=0)
            px_rate_list_background = torch.cat(px_rate_list_background, dim=0)
        if give_mean is True and n_samples > 1:
            px_rate_list_foreground = torch.mean(px_rate_list_foreground, dim=-1)
            px_rate_list_foreground = px_rate_list_foreground.cpu().numpy()
            px_rate_list_background = torch.mean(px_rate_list_background, dim=-1)
            px_rate_list_background = px_rate_list_background.cpu().numpy()
        return px_rate_list_foreground, px_rate_list_background

    @torch.no_grad()
    def imputation_list(self, n_samples: int = 1):
        original_list = []
        imputed_list = []
        batch_size = self.data_loader_kwargs["batch_size"] // n_samples
        for tensors, corrupted_tensors in zip(
            self.uncorrupted().sequential(batch_size=batch_size),
            self.corrupted().sequential(batch_size=batch_size),
        ):
            batch = tensors[0]
            actual_batch_size = batch.size(0)
            dropout_x, _, _, batch_index, labels, y = corrupted_tensors
            px_rate = self.model.get_sample_rate(
                dropout_x, y, batch_index=batch_index, label=labels, n_samples=n_samples
            )
            px_rate = px_rate[:, : self.gene_dataset.nb_genes]

            indices_dropout = torch.nonzero(batch - dropout_x)
            if indices_dropout.size() != torch.Size([0]):
                i = indices_dropout[:, 0]
                j = indices_dropout[:, 1]

                batch = batch.unsqueeze(0).expand(
                    (n_samples, batch.size(0), batch.size(1))
                )
                original = np.array(batch[:, i, j].view(-1).cpu())
                imputed = np.array(px_rate[..., i, j].view(-1).cpu())

                cells_index = np.tile(np.array(i.cpu()), n_samples)

                original_list += [
                    original[cells_index == i] for i in range(actual_batch_size)
                ]
                imputed_list += [
                    imputed[cells_index == i] for i in range(actual_batch_size)
                ]
            else:
                original_list = np.array([])
                imputed_list = np.array([])
        return original_list, imputed_list

    @torch.no_grad()
    def generate_parameters(self):
        raise NotImplementedError


class TotalTrainer(UnsupervisedTrainer):
    r"""The VariationalInference class for the unsupervised training of an autoencoder.

    Args:
        :model: A model instance from class ``TOTALVI``
        :gene_dataset: A gene_dataset instance like ``CbmcDataset()`` with attribute ``protein_expression``
        :train_size: The train size, either a float between 0 and 1 or and integer for the number of training samples
         to use Default: ``0.93``.
        :test_size: The test size, either a float between 0 and 1 or and integer for the number of training samples
         to use Default: ``0.02``. Note that if train and test do not add to 1 the remainder is placed in a validation set
        :\*\*kwargs: Other keywords arguments from the general Trainer class.
    """
    default_metrics_to_monitor = ["elbo"]

    def __init__(
        self,
        model,
        dataset,
        train_size=0.93,
        test_size=0.02,
        pro_weight=1.0,
        back_kl_weight=1.0,
        n_epochs_back_kl_warmup=200,
        n_epochs_kl_warmup=200,
        n_epochs_gene_warmup=None,
        **kwargs
    ):
        self.n_genes = dataset.nb_genes
        self.n_proteins = model.n_input_proteins

        self.pro_weight = pro_weight
        self.n_epochs_gene_warmup = n_epochs_gene_warmup
        self.n_epochs_back_kl_warmup = n_epochs_back_kl_warmup
        self.back_kl_weight = back_kl_weight
        super().__init__(
            model, dataset, n_epochs_kl_warmup=n_epochs_kl_warmup, **kwargs
        )
        if type(self) is TotalTrainer:
            self.train_set, self.test_set, self.validation_set = self.train_test_validation(
                model, dataset, train_size, test_size, type_class=TotalPosterior
            )
            self.train_set.to_monitor = []
            self.test_set.to_monitor = ["elbo"]
            self.validation_set.to_monitor = ["elbo"]

    def loss(self, tensors):
        sample_batch_X, local_l_mean, local_l_var, batch_index, label, sample_batch_Y = (
            tensors
        )
        reconst_loss_gene, reconst_loss_protein, kl_divergence, back_kl = self.model(
            sample_batch_X,
            sample_batch_Y,
            local_l_mean,
            local_l_var,
            batch_index,
            label,
        )

        loss = torch.mean(
            self.gene_weight * reconst_loss_gene
            + self.pro_weight * reconst_loss_protein
            + self.kl_weight * kl_divergence
            + self.back_kl_weight * self.back_warmup_weight * back_kl
        )
        return loss

    def on_epoch_begin(self):
        super().on_epoch_begin()
        if self.n_epochs_back_kl_warmup is not None:
            self.back_warmup_weight = min(1, self.epoch / self.n_epochs_back_kl_warmup)
        else:
            self.back_warmup_weight = 1.0
        if self.n_epochs_gene_warmup is not None:
            self.gene_weight = min(1, self.epoch / self.n_epochs_gene_warmup)
        else:
            self.gene_weight = 1.0
