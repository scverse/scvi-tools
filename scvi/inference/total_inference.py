from typing import Optional, Union, List, Callable, Tuple
import logging
import torch
from torch.distributions import Poisson, Gamma, Bernoulli, Normal
from torch.utils.data import DataLoader
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from scvi.inference import Posterior
from . import UnsupervisedTrainer

from scvi.dataset import GeneExpressionDataset
from scvi.models import TOTALVI, Classifier
from scvi.models.utils import one_hot

logger = logging.getLogger(__name__)


class TotalPosterior(Posterior):
    """The functional data unit for totalVI.

    A `TotalPosterior` instance is instantiated with a model and
    a `gene_dataset`, and as well as additional arguments that for Pytorch's `DataLoader`. A subset of indices
    can be specified, for purposes such as splitting the data into train/test/validation. Each trainer instance of the `TotalTrainer` class can therefore have multiple
    `TotalPosterior` instances to train a model. A `TotalPosterior` instance also comes with many methods or
    utilities for its corresponding data.


    :param model: A model instance from class ``TOTALVI``
    :param gene_dataset: A gene_dataset instance like ``CbmcDataset()`` with attribute ``protein_expression``
    :param shuffle: Specifies if a `RandomSampler` or a `SequentialSampler` should be used
    :param indices: Specifies how the data should be split with regards to train/test or labelled/unlabelled
    :param use_cuda: Default: ``True``
    :param data_loader_kwargs: Keyword arguments to passed into the `DataLoader`

    Examples:

    Let us instantiate a `trainer`, with a gene_dataset and a model

        >>> gene_dataset = CbmcDataset()
        >>> totalvi = TOTALVI(gene_dataset.nb_genes, len(gene_dataset.protein_names),
        ... n_batch=gene_dataset.n_batches, use_cuda=True)
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
        """Computes the ELBO.

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
                **kwargs,
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
                **kwargs,
            )
            log_lkl_gene += torch.sum(reconst_loss_gene).item()
            log_lkl_protein += torch.sum(reconst_loss_protein).item()

        n_samples = len(self.indices)
        return log_lkl_gene / n_samples, log_lkl_protein / n_samples

    def compute_marginal_log_likelihood(
        self, n_samples_mc: int = 100, batch_size: int = 96
    ):
        """Computes a biased estimator for log p(x, y), which is the marginal log likelihood.

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
    def get_latent(
        self, sample: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Output posterior z mean or sample, batch index, and label

        :param sample: z mean or z sample
        :return: 4-tuple of latent, batch_indices, labels, library_gene
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
    def generate(
        self, n_samples: int = 100, batch_size: int = 64
    ) -> Tuple[
        np.ndarray, np.ndarray
    ]:  # with n_samples>1 return original list/ otherwise sequential
        """Sample from posterior predictive. Proteins are concatenated to genes.

        :param n_samples: Number of posterior predictive samples
        :return: Tuple of posterior samples, original data
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
            px_dropout = torch.sigmoid(outputs["px_"]["dropout"])
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
    def get_sample_mixing(
        self,
        n_samples: int = 1,
        give_mean: bool = True,
        transform_batch: Optional[Union[int, List[int]]] = None,
    ) -> np.ndarray:
        """ Returns mixing bernoulli parameter for protein negative binomial mixtures (probability background)

        :param n_samples: number of samples from posterior distribution
        :param sample_protein_mixing: Sample mixing bernoulli, setting background to zero
        :param give_mean: bool, whether to return samples along first axis or average over samples
        :param transform_batch: Batches to condition on.
        If transform_batch is:
            - None, then real observed batch is used
            - int, then batch transform_batch is used
            - list of int, then values are averaged over provided batches.
        :return: array of probability background
        """
        py_mixings = []
        if (transform_batch is None) or (isinstance(transform_batch, int)):
            transform_batch = [transform_batch]
        for tensors in self:
            x, _, _, batch_index, label, y = tensors
            py_mixing = torch.zeros_like(y)
            if n_samples > 1:
                py_mixing = torch.stack(n_samples * [py_mixing])
            for b in transform_batch:
                outputs = self.model.inference(
                    x,
                    y,
                    batch_index=batch_index,
                    label=label,
                    n_samples=n_samples,
                    transform_batch=b,
                )
                py_mixing += torch.sigmoid(outputs["py_"]["mixing"])
            py_mixing /= len(transform_batch)
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
    def get_sample_scale(
        self,
        transform_batch=None,
        eps=0.5,
        normalize_pro=False,
        sample_bern=True,
        include_bg=False,
    ) -> np.ndarray:
        """Helper function to provide normalized expression for DE testing.

        For normalized, denoised expression, please use
            `get_normalized_denoised_expression()`

        :param transform_batch: Int of batch to "transform" all cells into
        :param eps: Prior count to add to protein normalized expression
        :param normalize_pro: bool, whether to make protein expression sum to one in a cell
        :param include_bg: bool, whether to include the background component of expression
        """
        scales = []
        for tensors in self:
            x, _, _, batch_index, label, y = tensors
            model_scale = self.model.get_sample_scale(
                x,
                y,
                batch_index=batch_index,
                label=label,
                n_samples=1,
                transform_batch=transform_batch,
                eps=eps,
                normalize_pro=normalize_pro,
                sample_bern=sample_bern,
                include_bg=include_bg,
            )
            # prior count for proteins
            scales += [torch.cat(model_scale, dim=-1).cpu().numpy()]
        return np.concatenate(scales)

    @torch.no_grad()
    def get_normalized_denoised_expression(
        self,
        n_samples: int = 1,
        give_mean: bool = True,
        transform_batch: Optional[Union[int, List[int]]] = None,
        sample_protein_mixing: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Returns the tensors of denoised normalized gene and protein expression

        :param n_samples: number of samples from posterior distribution
        :param sample_protein_mixing: Sample mixing bernoulli, setting background to zero
        :param give_mean: bool, whether to return samples along first axis or average over samples
        :param transform_batch: Batches to condition on.
        If transform_batch is:
            - None, then real observed batch is used
            - int, then batch transform_batch is used
            - list of int, then values are averaged over provided batches.
        :return: Denoised genes, denoised proteins
        """

        scale_list_gene = []
        scale_list_pro = []
        if (transform_batch is None) or (isinstance(transform_batch, int)):
            transform_batch = [transform_batch]
        for tensors in self:
            x, _, _, batch_index, label, y = tensors
            px_scale = torch.zeros_like(x)
            py_scale = torch.zeros_like(y)
            if n_samples > 1:
                px_scale = torch.stack(n_samples * [px_scale])
                py_scale = torch.stack(n_samples * [py_scale])
            for b in transform_batch:
                outputs = self.model.inference(
                    x,
                    y,
                    batch_index=batch_index,
                    label=label,
                    n_samples=n_samples,
                    transform_batch=b,
                )
                px_scale += outputs["px_"]["scale"]

                py_ = outputs["py_"]
                # probability of background
                protein_mixing = 1 / (1 + torch.exp(-py_["mixing"]))
                if sample_protein_mixing is True:
                    protein_mixing = Bernoulli(protein_mixing).sample()
                py_scale += py_["rate_fore"] * (1 - protein_mixing)
            px_scale /= len(transform_batch)
            py_scale /= len(transform_batch)
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
    def get_protein_mean(
        self,
        n_samples: int = 1,
        give_mean: bool = True,
        transform_batch: Optional[Union[int, List[int]]] = None,
    ) -> np.ndarray:
        """Returns the tensors of protein mean (with foreground and background)

        :param n_samples: number of samples from posterior distribution
        :param give_mean: bool, whether to return samples along first axis or average over samples
        :param transform_batch: Batches to condition on.
        If transform_batch is:
            - None, then real observed batch is used
            - int, then batch transform_batch is used
            - list of int, then values are averaged over provided batches.
        """
        if (transform_batch is None) or (isinstance(transform_batch, int)):
            transform_batch = [transform_batch]
        rate_list_pro = []
        for tensors in self:
            x, _, _, batch_index, label, y = tensors
            protein_rate = torch.zeros_like(y)
            if n_samples > 1:
                protein_rate = torch.stack(n_samples * [protein_rate])
            for b in transform_batch:
                outputs = self.model.inference(
                    x,
                    y,
                    batch_index=batch_index,
                    label=label,
                    n_samples=n_samples,
                    transform_batch=b,
                )
                py_ = outputs["py_"]
                pi = 1 / (1 + torch.exp(-py_["mixing"]))
                protein_rate += py_["rate_fore"] * (1 - pi) + py_["rate_back"] * pi
            protein_rate /= len(transform_batch)
            rate_list_pro.append(protein_rate.cpu())

        if n_samples > 1:
            # concatenate along batch dimension -> result shape = (samples, cells, features)
            rate_list_pro = torch.cat(rate_list_pro, dim=1)
            # (cells, features, samples)
            rate_list_pro = rate_list_pro.permute(1, 2, 0)
        else:
            rate_list_pro = torch.cat(rate_list_pro, dim=0)

        if give_mean is True and n_samples > 1:
            rate_list_pro = torch.mean(rate_list_pro, dim=-1)

        rate_list_pro = rate_list_pro.cpu().numpy()

        return rate_list_pro

    @torch.no_grad()
    def generate_denoised_samples(
        self,
        n_samples: int = 25,
        batch_size: int = 64,
        rna_size_factor: int = 1,
        transform_batch: Optional[int] = None,
    ):
        """Samples from an adjusted posterior predictive. Proteins are concatenated to genes.

        :param n_samples: How may samples per cell
        :param batch_size: Mini-batch size for sampling. Lower means less GPU memory footprint
        :rna_size_factor: size factor for RNA prior to sampling gamma distribution
        :transform_batch: int of which batch to condition on for all cells
        :return:
        """
        posterior_list = []
        for tensors in self.update({"batch_size": batch_size}):
            x, _, _, batch_index, labels, y = tensors
            with torch.no_grad():
                outputs = self.model.inference(
                    x,
                    y,
                    batch_index=batch_index,
                    label=labels,
                    n_samples=n_samples,
                    transform_batch=transform_batch,
                )
            px_ = outputs["px_"]
            py_ = outputs["py_"]

            pi = 1 / (1 + torch.exp(-py_["mixing"]))
            mixing_sample = Bernoulli(pi).sample()
            protein_rate = py_["rate_fore"]
            rate = torch.cat((rna_size_factor * px_["scale"], protein_rate), dim=-1)
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
            data = l_train.cpu().numpy()
            # make background 0
            data[:, :, self.gene_dataset.nb_genes :] = (
                data[:, :, self.gene_dataset.nb_genes :]
                * (1 - mixing_sample).cpu().numpy()
            )
            posterior_list += [data]

            posterior_list[-1] = np.transpose(posterior_list[-1], (1, 2, 0))

        return np.concatenate(posterior_list, axis=0)

    @torch.no_grad()
    def generate_feature_correlation_matrix(
        self,
        n_samples: int = 25,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[Union[int, List[int]]] = None,
        correlation_mode: str = "pearson",
        log_transform: bool = False,
    ):
        """ Wrapper of `generate_denoised_samples()` to create a gene-protein gene-protein corr matrix
        :param n_samples: How may samples per cell
        :param batch_size: Mini-batch size for sampling. Lower means less GPU memory footprint
        :rna_size_factor: size factor for RNA prior to sampling gamma distribution
        :param transform_batch: Batches to condition on.
        If transform_batch is:
            - None, then real observed batch is used
            - int, then batch transform_batch is used
            - list of int, then values are averaged over provided batches.
        :param log_transform: Whether to log transform denoised values prior to correlation calculation
        :return: Correlation matrix
        """
        if (transform_batch is None) or (isinstance(transform_batch, int)):
            transform_batch = [transform_batch]
        corr_mats = []
        for b in transform_batch:
            denoised_data = self.generate_denoised_samples(
                n_samples=n_samples,
                batch_size=batch_size,
                rna_size_factor=rna_size_factor,
                transform_batch=b,
            )
            flattened = np.zeros(
                (denoised_data.shape[0] * n_samples, denoised_data.shape[1])
            )
            for i in range(n_samples):
                flattened[
                    denoised_data.shape[0] * (i) : denoised_data.shape[0] * (i + 1)
                ] = denoised_data[:, :, i]
            if log_transform is True:
                flattened[:, : self.gene_dataset.nb_genes] = np.log(
                    flattened[:, : self.gene_dataset.nb_genes] + 1e-8
                )
                flattened[:, self.gene_dataset.nb_genes :] = np.log1p(
                    flattened[:, self.gene_dataset.nb_genes :]
                )
            if correlation_mode == "pearson":
                corr_matrix = np.corrcoef(flattened, rowvar=False)
            else:
                corr_matrix = spearmanr(flattened, axis=0)
            corr_mats.append(corr_matrix)
        corr_matrix = np.mean(np.stack(corr_mats), axis=0)
        return corr_matrix

    @torch.no_grad()
    def imputation(self, n_samples: int = 1):
        """Gene imputation
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
        """This code is identical to same function in posterior.py

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
    def differential_expression_score(
        self,
        idx1: Union[List[bool], np.ndarray],
        idx2: Union[List[bool], np.ndarray],
        mode: Optional[str] = "vanilla",
        batchid1: Optional[Union[List[int], np.ndarray]] = None,
        batchid2: Optional[Union[List[int], np.ndarray]] = None,
        use_observed_batches: Optional[bool] = False,
        n_samples: int = 5000,
        use_permutation: bool = True,
        M_permutation: int = 10000,
        all_stats: bool = True,
        change_fn: Optional[Union[str, Callable]] = None,
        m1_domain_fn: Optional[Callable] = None,
        delta: Optional[float] = 0.5,
        cred_interval_lvls: Optional[Union[List[float], np.ndarray]] = None,
        **kwargs,
    ) -> pd.DataFrame:
        r"""Unified method for differential expression inference.

        This function is an extension of the `get_bayes_factors` method
        providing additional genes information to the user

        Two modes coexist:

        - the "vanilla" mode follows protocol described in [Lopez18]_
        In this case, we perform hypothesis testing based on the hypotheses

        .. math::
            M_1: h_1 > h_2 ~\text{and}~ M_2: h_1 \leq h_2

        DE can then be based on the study of the Bayes factors

        .. math::
            \log p(M_1 | x_1, x_2) / p(M_2 | x_1, x_2)

        - the "change" mode (described in [Boyeau19]_)
        consists in estimating an effect size random variable (e.g., log fold-change) and
        performing Bayesian hypothesis testing on this variable.
        The `change_fn` function computes the effect size variable r based two inputs
        corresponding to the normalized means in both populations.

        Hypotheses:

        .. math::
            M_1: r \in R_1 ~\text{(effect size r in region inducing differential expression)}

        .. math::
            M_2: r  \notin R_1 ~\text{(no differential expression)}

        To characterize the region :math:`R_1`, which induces DE, the user has two choices.

        1. A common case is when the region :math:`[-\delta, \delta]` does not induce differential
        expression.
        If the user specifies a threshold delta,
        we suppose that :math:`R_1 = \mathbb{R} \setminus [-\delta, \delta]`

        2. specify an specific indicator function

        .. math::
            f: \mathbb{R} \mapsto \{0, 1\} ~\text{s.t.}~ r \in R_1 ~\text{iff.}~ f(r) = 1

        Decision-making can then be based on the estimates of

        .. math::
            p(M_1 \mid x_1, x_2)

        Both modes require to sample the normalized means posteriors.
        To that purpose, we sample the Posterior in the following way:

        1. The posterior is sampled n_samples times for each subpopulation

        2. For computation efficiency (posterior sampling is quite expensive), instead of
            comparing the obtained samples element-wise, we can permute posterior samples.
            Remember that computing the Bayes Factor requires sampling
            :math:`q(z_A \mid x_A)` and :math:`q(z_B \mid x_B)`

        Currently, the code covers several batch handling configurations:

        1. If ``use_observed_batches=True``, then batch are considered as observations
        and cells' normalized means are conditioned on real batch observations

        2. If case (cell group 1) and control (cell group 2) are conditioned on the same
        batch ids.
        Examples:
            >>> set(batchid1) = set(batchid2)

        or
            >>> batchid1 = batchid2 = None


        3. If case and control are conditioned on different batch ids that do not intersect
        i.e.,
            >>> set(batchid1) != set(batchid2)

        and
            >>> len(set(batchid1).intersection(set(batchid2))) == 0

        This function does not cover other cases yet and will warn users in such cases.

        :param mode: one of ["vanilla", "change"]

        :param idx1: bool array masking subpopulation cells 1. Should be True where cell is
          from associated population
        :param idx2: bool array masking subpopulation cells 2. Should be True where cell is
          from associated population
        :param batchid1: List of batch ids for which you want to perform DE Analysis for
          subpopulation 1. By default, all ids are taken into account
        :param batchid2: List of batch ids for which you want to perform DE Analysis for
          subpopulation 2. By default, all ids are taken into account
        :param use_observed_batches: Whether normalized means are conditioned on observed
          batches

        :param n_samples: Number of posterior samples
        :param use_permutation: Activates step 2 described above.
          Simply formulated, pairs obtained from posterior sampling (when calling
          `sample_scale_from_batch`) will be randomly permuted so that the number of
          pairs used to compute Bayes Factors becomes M_permutation.
        :param M_permutation: Number of times we will "mix" posterior samples in step 2.
          Only makes sense when use_permutation=True

        :param change_fn: function computing effect size based on both normalized means
        :param m1_domain_fn: custom indicator function of effect size regions
          inducing differential expression
        :param delta: specific case of region inducing differential expression.
          In this case, we suppose that :math:`R \setminus [-\delta, \delta]` does not induce differential expression
          (LFC case)
        :param cred_interval_lvls: List of credible interval levels to compute for the posterior
          LFC distribution

        :param all_stats: whether additional metrics should be provided
        :\**kwargs: Other keywords arguments for `get_sample_scale()`

        :return: Differential expression properties
        """
        all_info = self.get_bayes_factors(
            idx1=idx1,
            idx2=idx2,
            mode=mode,
            batchid1=batchid1,
            batchid2=batchid2,
            use_observed_batches=use_observed_batches,
            n_samples=n_samples,
            use_permutation=use_permutation,
            M_permutation=M_permutation,
            change_fn=change_fn,
            m1_domain_fn=m1_domain_fn,
            cred_interval_lvls=cred_interval_lvls,
            delta=delta,
            **kwargs,
        )
        col_names = np.concatenate(
            [self.gene_dataset.gene_names, self.gene_dataset.protein_names]
        )
        if all_stats is True:
            nan = np.array([np.nan] * len(self.gene_dataset.protein_names))
            (
                mean1,
                mean2,
                nonz1,
                nonz2,
                norm_mean1,
                norm_mean2,
            ) = self.gene_dataset.raw_counts_properties(idx1, idx2)
            mean1_pro = self.gene_dataset.protein_expression[idx1, :].mean(0)
            mean2_pro = self.gene_dataset.protein_expression[idx2, :].mean(0)
            nonz1_pro = (self.gene_dataset.protein_expression[idx1, :] > 0).mean(0)
            nonz2_pro = (self.gene_dataset.protein_expression[idx2, :] > 0).mean(0)
            # TODO implement properties for proteins
            genes_properties_dict = dict(
                raw_mean1=np.concatenate([mean1, mean1_pro]),
                raw_mean2=np.concatenate([mean2, mean2_pro]),
                non_zeros_proportion1=np.concatenate([nonz1, nonz1_pro]),
                non_zeros_proportion2=np.concatenate([nonz2, nonz2_pro]),
                raw_normalized_mean1=np.concatenate([norm_mean1, nan]),
                raw_normalized_mean2=np.concatenate([norm_mean2, nan]),
            )
            all_info = {**all_info, **genes_properties_dict}

        res = pd.DataFrame(all_info, index=col_names)
        sort_key = "proba_de" if mode == "change" else "bayes_factor"
        res = res.sort_values(by=sort_key, ascending=False)
        return res

    @torch.no_grad()
    def generate_parameters(self):
        raise NotImplementedError


default_early_stopping_kwargs = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "patience": 45,
    "threshold": 0,
    "reduce_lr_on_plateau": True,
    "lr_patience": 30,
    "lr_factor": 0.6,
    "posterior_class": TotalPosterior,
}


class TotalTrainer(UnsupervisedTrainer):
    r"""Unsupervised training for totalVI using variational inference

    :param model: A model instance from class ``TOTALVI``
    :param gene_dataset: A gene_dataset instance like ``CbmcDataset()`` with attribute ``protein_expression``
    :param train_size: The train size, either a float between 0 and 1 or and integer for the number of training samples
        to use Default: ``0.90``.
    :param test_size: The test size, either a float between 0 and 1 or and integer for the number of training samples
        to use Default: ``0.10``. Note that if train and test do not add to 1 the remainder is placed in a validation set
    :param pro_recons_weight: Scaling factor on the reconstruction loss for proteins. Default: ``1.0``.
    :param n_epochs_kl_warmup: Number of epochs for annealing the KL terms for `z` and `mu` of the ELBO (from 0 to 1). If None, no warmup performed, unless
        `n_iter_kl_warmup` is set.
    :param n_iter_kl_warmup: Number of minibatches for annealing the KL terms for `z` and `mu` of the ELBO (from 0 to 1). If set to "auto", the number
        of iterations is equal to 75% of the number of cells. `n_epochs_kl_warmup` takes precedence if it is not None. If both are None, then
        no warmup is performed.
    :param discriminator: Classifier used for adversarial training scheme
    :param use_adversarial_loss: Whether to use adversarial classifier to improve mixing
    :param kappa: Scaling factor for adversarial loss. If None, follow inverse of kl warmup schedule.
    :param early_stopping_kwargs: Keyword args for early stopping. If "auto", use totalVI defaults. If None, disable early stopping.
    :param \*\*kwargs: Other keywords arguments from the general Trainer class.
    """

    default_metrics_to_monitor = ["elbo"]

    def __init__(
        self,
        model: TOTALVI,
        dataset: GeneExpressionDataset,
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
        self.n_genes = dataset.nb_genes
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
                n_labels=self.gene_dataset.n_batches,
                n_layers=2,
                logits=True,
            )

        self.discriminator = discriminator
        if self.use_cuda and self.discriminator is not None:
            self.discriminator.cuda()

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
            + self.kl_weight * kl_div_back_pro
        )

        return loss

    def loss_discriminator(
        self, z, batch_index, predict_true_class=True, return_details=True
    ):

        n_classes = self.gene_dataset.n_batches
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
            sample_batch_X,
            local_l_mean,
            local_l_var,
            batch_index,
            label,
            sample_batch_Y,
        ) = tensors

        z = self.model.sample_from_posterior_z(
            sample_batch_X, sample_batch_Y, batch_index, give_mean=False
        )

        return z

    def train(self, n_epochs=500, lr=4e-3, eps=0.01, params=None):

        super().train(n_epochs=n_epochs, lr=lr, eps=eps, params=params)

    def on_training_loop(self, tensors_list):
        if self.use_adversarial_loss:
            if self.kappa is None:
                kappa = 1 - self.kl_weight
            else:
                kappa = self.kappa
            batch_index = tensors_list[0][3]
            if kappa > 0:
                z = self._get_z(*tensors_list)
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
            self.current_loss = loss = self.loss(*tensors_list)
            if kappa > 0:
                (loss + fool_loss).backward()
            else:
                loss.backward()
            self.optimizer.step()

        else:
            self.current_loss = loss = self.loss(*tensors_list)
            self.optimizer.zero_grad()
            loss.backward()
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
