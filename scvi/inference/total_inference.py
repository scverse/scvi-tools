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
    can be specified, for purposes such as splitting the data into train/test/validation. Each trainer instance of the `TotalTrainer` class can therefore have multiple
    `TotalPosterior` instances to train a model. A `TotalPosterior` instance also comes with many methods or
    utilities for its corresponding data.


    :param model: A model instance from class ``TOTALVI``
    :param gene_dataset: A gene_dataset instance like ``CbmcDataset()`` with attribute ``protein_expression``
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
    def reconstruction_error(self, mode="total"):
        ll_gene, ll_protein = self.compute_reconstruction_error(self.model)
        if mode == "total":
            return ll_gene + ll_protein
        elif mode == "gene":
            return ll_gene
        else:
            return ll_protein

    reconstruction_error.mode = "min"

    @torch.no_grad()
    def marginal_ll(self, n_mc_samples=1000):
        ll = self.compute_marginal_log_likelihood()
        return ll

    @torch.no_grad()
    def get_protein_background_mean(self):
        background_mean = []
        for tensors in self:
            x, _, _, batch_index, label, y = tensors
            outputs = self.model.inference(
                x, y, batch_index=batch_index, label=label, n_samples=1
            )
            b_mean = outputs["py_"]["rate_back"]
            background_mean += [np.array(b_mean.cpu())]
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
            (
                reconst_loss_gene,
                reconst_loss_protein,
                kl_div_z,
                kl_div_gene_l,
                kl_div_back_pro,
            ) = vae(
                x,
                y,
                local_l_mean,
                local_l_var,
                batch_index=batch_index,
                label=labels,
                **kwargs
            )
            elbo += torch.sum(
                reconst_loss_gene
                + reconst_loss_protein
                + kl_div_z
                + kl_div_gene_l
                + kl_div_back_pro
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
            (
                reconst_loss_gene,
                reconst_loss_protein,
                kl_div_z,
                kl_div_l_gene,
                kl_div_back_pro,
            ) = vae(
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
        (a fairly high value like 100 should be enough). 5000 is the standard in machine learning publications.
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
                px_ = outputs["px_"]
                py_ = outputs["py_"]
                log_library = outputs["untran_l"]
                # really need not softmax transformed random variable
                z = outputs["untran_z"]
                log_pro_back_mean = outputs["log_pro_back_mean"]

                # Reconstruction Loss
                (
                    reconst_loss_gene,
                    reconst_loss_protein,
                ) = self.model.get_reconstruction_loss(x, y, px_, py_)

                # Log-probabilities
                p_l_gene = (
                    Normal(local_l_mean, local_l_var.sqrt())
                    .log_prob(log_library)
                    .sum(dim=-1)
                )
                p_z = Normal(0, 1).log_prob(z).sum(dim=-1)
                p_mu_back = self.model.back_mean_prior.log_prob(log_pro_back_mean).sum(
                    dim=-1
                )
                p_xy_zl = -(reconst_loss_gene + reconst_loss_protein)
                q_z_x = Normal(qz_m, qz_v.sqrt()).log_prob(z).sum(dim=-1)
                q_l_x = Normal(ql_m, ql_v.sqrt()).log_prob(log_library).sum(dim=-1)
                q_mu_back = (
                    Normal(py_["back_alpha"], py_["back_beta"])
                    .log_prob(log_pro_back_mean)
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
        :return: 4-tuple of np.ndarrays, latent, batch_indices, labels, library_gene
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
        scales = []
        fixed_batch = float(fixed_batch)
        for tensors in self:
            x, local_l_mean, local_l_var, batch_index, label, y = tensors
            scales += [
                torch.cat(self.model.scale_from_z(x, y, fixed_batch).cpu(), dim=-1)
            ]
        return np.concatenate(scales)

    @torch.no_grad()
    def generate(
        self,
        n_samples: int = 100,
        genes: Optional[np.ndarray] = None,
        batch_size: int = 64,
    ):  # with n_samples>1 return original list/ otherwise sequential
        """
        Return samples from posterior predictive. Proteins are concatenated to genes.
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
            px_ = outputs["px_"]
            py_ = outputs["py_"]

            pi = 1 / (1 + torch.exp(-py_["mixing"]))
            mixing_sample = Bernoulli(pi).sample()
            protein_rate = (
                py_["rate_fore"] * (1 - mixing_sample)
                + py_["rate_back"] * mixing_sample
            )
            rate = torch.cat((px_["rate"], protein_rate), dim=-1)
            if len(px_["r"].size()) == 2:
                px_dispersion = px_["r"]
            else:
                px_dispersion = torch.ones_like(x) * px_["r"]
            if len(py_["r"].size()) == 2:
                py_dispersion = py_["r"]
            else:
                py_dispersion = torch.ones_like(y) * py_["r"]

            dispersion = torch.cat((px_dispersion, py_dispersion), dim=-1)

            # This gamma is really l*w using scVI manuscript notation
            p = rate / (rate + dispersion)
            r = dispersion
            l_train = Gamma(r, (1 - p) / p).sample()
            data = Poisson(l_train).sample().cpu().numpy()
            # """
            # In numpy (shape, scale) => (concentration, rate), with scale = p /(1 - p)
            # rate = (1 - p) / p  # = 1/scale # used in pytorch
            # """
            original_list += [np.array(torch.cat((x, y), dim=-1).cpu())]
            posterior_list += [data]

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
        """ Zero-inflation mixing component for genes
        """
        px_dropouts = []
        for tensors in self:
            x, _, _, batch_index, label, y = tensors
            outputs = self.model.inference(
                x, y, batch_index=batch_index, label=label, n_samples=n_samples
            )
            px_dropout = outputs["px_"]["dropout"]
            px_dropouts += [px_dropout.cpu()]
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
    def get_sample_mixing(self, n_samples: int = 1, give_mean: bool = True):
        """ Returns mixing bernoulli parameter for negative binomial mixtures
        """
        py_mixings = []
        for tensors in self:
            x, _, _, batch_index, label, y = tensors
            outputs = self.model.inference(
                x, y, batch_index=batch_index, label=label, n_samples=n_samples
            )
            py_mixing = outputs["py_"]["mixing"]
            py_mixings += [py_mixing.cpu()]
        if n_samples > 1:
            # concatenate along batch dimension -> result shape = (samples, cells, features)
            py_mixings = torch.cat(py_mixings, dim=1)
            # (cells, features, samples)
            py_mixings = py_mixings.permute(1, 2, 0)
        else:
            py_mixings = torch.cat(py_mixings, dim=0)

        if give_mean is True and n_samples > 1:
            py_mixings = torch.mean(py_mixings, dim=-1)

        py_mixings = py_mixings.cpu().numpy()

        return py_mixings

    @torch.no_grad()
    def get_normalized_denoised_expression(
        self, n_samples: int = 1, give_mean: bool = True
    ):
        """Returns the tensors of denoised normalized gene and protein expression

        :param n_samples: number of samples from posterior distribution
        :param give_mean: bool, whether to return samples along first axis or average over samples
        :rtype: 2-tuple of :py:class:`np.ndarray`
        """

        scale_list_gene = []
        scale_list_pro = []
        for tensors in self:
            x, _, _, batch_index, label, y = tensors
            outputs = self.model.inference(
                x, y, batch_index=batch_index, label=label, n_samples=n_samples
            )
            px_scale = outputs["px_"]["scale"]
            py_scale = outputs["py_"]["scale"]

            scale_list_gene.append(px_scale.cpu())
            scale_list_pro.append(py_scale.cpu())

        if n_samples > 1:
            # concatenate along batch dimension -> result shape = (samples, cells, features)
            scale_list_gene = torch.cat(scale_list_gene, dim=1)
            scale_list_pro = torch.cat(scale_list_pro, dim=1)
            # (cells, features, samples)
            scale_list_gene = scale_list_gene.permute(1, 2, 0)
            scale_list_pro = scale_list_pro.permute(1, 2, 0)
        else:
            scale_list_gene = torch.cat(scale_list_gene, dim=0)
            scale_list_pro = torch.cat(scale_list_pro, dim=0)

        if give_mean is True and n_samples > 1:
            scale_list_gene = torch.mean(scale_list_gene, dim=-1)
            scale_list_pro = torch.mean(scale_list_pro, dim=-1)

        scale_list_gene = scale_list_gene.cpu().numpy()
        scale_list_pro = scale_list_pro.cpu().numpy()

        return scale_list_gene, scale_list_pro

    @torch.no_grad()
    def imputation(self, n_samples: int = 1):
        """ Gene imputation
        """
        imputed_list = []
        for tensors in self:
            x, _, _, batch_index, label, y = tensors
            px_rate = self.model.get_sample_rate(
                x, y, batch_index=batch_index, label=label, n_samples=n_samples
            )
            imputed_list += [np.array(px_rate.cpu())]
        imputed_list = np.concatenate(imputed_list)
        return imputed_list.squeeze()

    @torch.no_grad()
    def imputation_list(self, n_samples: int = 1):
        """ This code is identical to same function in posterior.py

            Except, we use the totalVI definition of `model.get_sample_rate`
        """
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

    @torch.no_grad()
    def get_sample_scale(self):
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
        train_size=0.90,
        test_size=0.05,
        pro_recons_weight=1.0,
        n_epochs_back_kl_warmup=200,
        n_epochs_kl_warmup=200,
        **kwargs
    ):
        self.n_genes = dataset.nb_genes
        self.n_proteins = model.n_input_proteins

        self.pro_recons_weight = pro_recons_weight
        self.n_epochs_back_kl_warmup = n_epochs_back_kl_warmup
        super().__init__(
            model, dataset, n_epochs_kl_warmup=n_epochs_kl_warmup, **kwargs
        )
        if type(self) is TotalTrainer:
            (
                self.train_set,
                self.test_set,
                self.validation_set,
            ) = self.train_test_validation(
                model, dataset, train_size, test_size, type_class=TotalPosterior
            )
            self.train_set.to_monitor = []
            self.test_set.to_monitor = ["elbo"]
            self.validation_set.to_monitor = ["elbo"]

    def loss(self, tensors):
        (
            sample_batch_X,
            local_l_mean,
            local_l_var,
            batch_index,
            label,
            sample_batch_Y,
        ) = tensors
        (
            reconst_loss_gene,
            reconst_loss_protein,
            kl_div_z,
            kl_div_l_gene,
            kl_div_back_pro,
        ) = self.model(
            sample_batch_X,
            sample_batch_Y,
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
            + self.back_warmup_weight * kl_div_back_pro
        )
        return loss

    def on_epoch_begin(self):
        super().on_epoch_begin()
        if self.n_epochs_back_kl_warmup is not None:
            self.back_warmup_weight = min(1, self.epoch / self.n_epochs_back_kl_warmup)
        else:
            self.back_warmup_weight = 1.0
