import numpy as np
import logging
import torch
import pandas as pd
from anndata import AnnData

from typing import Optional, Union, List, Dict, Tuple
from scvi._compat import Literal
from scvi.models._modules.totalvae import TOTALVAE
from scvi.models import SCVI

# from scvi.models._differential import DifferentialExpression
from scvi import _CONSTANTS
from scvi.inference.total_inference import TotalPosterior, TotalTrainer

logger = logging.getLogger(__name__)


class TOTALVI(SCVI):
    """total Variational Inference [GayosoSteier20]_

    Parameters
    ----------
    adata
        AnnData object that has been registered with scvi
    n_latent
        Dimensionality of the latent space
    gene_dispersion
        One of the following

        * ``'gene'`` - genes_dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - genes_dispersion can differ between different batches
        * ``'gene-label'`` - genes_dispersion can differ between different labels
    protein_dispersion
        One of the following

        * ``'protein'`` - protein_dispersion parameter is constant per protein across cells
        * ``'protein-batch'`` - protein_dispersion can differ between different batches NOT TESTED
        * ``'protein-label'`` - protein_dispersion can differ between different labels NOT TESTED
    gene_likelihood
        One of

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
    latent_distribution
        One of

        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)

    Examples
    --------

    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.dataset.setup_anndata(adata, batch_key="batch", protein_expression_obsm_key="protein_expression")
    >>> vae = scvi.models.TOTALVI(adata)
    >>> vae.train(n_epochs=400)
    >>> adata.obsm["X_scVI"] = vae.get_latent_representation()
    """

    def __init__(
        self,
        adata: AnnData,
        n_latent: int = 20,
        gene_dispersion: Literal[
            "gene", "gene-batch", "gene-label", "gene-cell"
        ] = "gene",
        protein_dispersion: Literal[
            "protein", "protein-batch", "protein-label"
        ] = "protein",
        gene_likelihood: Literal["zinb", "nb"] = "nb",
        latent_distribution: Literal["normal", "ln"] = "ln",
        use_cuda: bool = True,
        **model_kwargs,
    ):
        assert (
            "scvi_data_registry" in adata.uns.keys()
        ), "Please setup your AnnData with scvi.dataset.setup_anndata(adata) first"

        self.adata = adata
        summary_stats = adata.uns["scvi_summary_stats"]

        if "totalvi_batch_mask" in adata.uns.keys():
            batch_mask = adata.uns["totalvi_batch_mask"]
        else:
            batch_mask = None
        self.model = TOTALVAE(
            n_input=summary_stats["n_genes"],
            n_input_proteins=summary_stats["n_proteins"],
            n_batch=summary_stats["n_batch"],
            gene_dispersion=gene_dispersion,
            protein_dispersion=protein_dispersion,
            reconstruction_loss_gene=gene_likelihood,
            latent_distribution=latent_distribution,
            protein_batch_mask=batch_mask,
            **model_kwargs,
        )
        self.is_trained = False
        self.use_cuda = use_cuda and torch.cuda.is_available()
        self.batch_size = 256
        self._posterior_class = TotalPosterior

    def train(
        self,
        n_epochs=400,
        train_size=0.9,
        test_size=None,
        lr=1e-3,
        n_iter_kl_warmup=None,
        n_epochs_kl_warmup=400,
        metric_frequency=1,
        trainer_kwargs={},
        train_kwargs={},
    ):

        if "totalvi_batch_mask" in self.adata.uns.keys():
            imputation = True
        self.trainer = TotalTrainer(
            self.model,
            self.adata,
            train_size=train_size,
            test_size=test_size,
            n_iter_kl_warmup=n_iter_kl_warmup,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            frequency=metric_frequency,
            batch_size=self.batch_size,
            use_adversarial_loss=imputation,
            **trainer_kwargs,
        )
        self.trainer.train(n_epochs=n_epochs, lr=lr, **train_kwargs)
        self.is_trained = True
        self.train_indices = self.trainer.train_set.indices
        self.test_indices = self.trainer.test_set.indices
        self.validation_indices = self.trainer.validation_set.indices

    @torch.no_grad()
    def get_reconstruction_error(self, adata=None, indices=None, mode="total"):

        post = self._make_posterior(adata=adata, indices=indices)

        return -post.reconstruction_error(mode=mode)

    @torch.no_grad()
    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Union[np.ndarray, List[int]]] = None,
        give_mean: bool = True,
        mc_samples: int = 5000,
    ) -> np.ndarray:
        """Return the latent representation for each cell

        Parameters
        ----------
        adata
            AnnData object that has been registered with scvi. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Give mean of distribution or sample from it
        mc_samples
            For distributions with no closed-form mean (e.g., `logistic normal`), how many Monte Carlo
            samples to take for computing mean.

        Returns
        -------
        latent_representation : np.ndarray
            Low-dimensional representation for each cell

        Examples
        --------

        >>> vae = scvi.model.TOTALVI(adata)
        >>> vae.train(n_epochs=400)
        >>> adata.obsm["X_totalVI"] = vae.get_latent_representation()

        We can also get the latent representation for a subset of cells

        >>> adata_subset = adata[adata.obs.cell_type == "really cool cell type"]
        >>> latent_subset = vae.get_latent_representation(adata_subset)

        """

        if self.is_trained is False:
            raise RuntimeError("Please train the model first.")

        post = self._make_posterior(adata=adata, indices=indices)
        latent = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
            batch = tensors[_CONSTANTS.BATCH_KEY]
            z = self.model.sample_from_posterior_z(
                x, y, batch, give_mean=give_mean, n_samples=mc_samples
            )
            latent += [z.cpu()]
        return np.array(torch.cat(latent))

    @torch.no_grad()
    def get_latent_library_size(self, adata=None, indices=None, give_mean=True):
        if self.is_trained is False:
            raise RuntimeError("Please train the model first.")

        post = self._make_posterior(adata=adata, indices=indices)
        libraries = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
            batch = tensors[_CONSTANTS.BATCH_KEY]
            library = self.model.sample_from_posterior_l(
                x, y, batch, give_mean=give_mean
            )
            libraries += [library.cpu()]
        return np.array(torch.cat(libraries))

    @torch.no_grad()
    def get_normalized_expression(
        self,
        adata=None,
        indices=None,
        transform_batch: Optional[int] = None,
        gene_list: Optional[Union[np.ndarray, List[int]]] = None,
        protein_list: Optional[Union[np.ndarray, List[int]]] = None,
        library_size: Optional[Union[float, Literal["latent"]]] = 1,
        n_samples: int = 1,
        sample_protein_mixing: bool = True,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
    ) -> Tuple[Union[np.ndarray, pd.DataFrame], Union[np.ndarray, pd.DataFrame]]:
        r"""Returns the normalized gene expression and protein expression

        This is denoted as :math:`\rho_n` in the totalVI paper for genes, and TODO
        for proteins, :math:`(1-\pi_{nt})\alpha_{nt}\beta_{nt}`.

        Parameters
        ----------
        transform_batch
            Batch to condition on.
            If transform_batch is:

            - None, then real observed batch is used
            - int, then batch transform_batch is used
            - List[int], then average over batches in list
        gene_list
            Return frequencies of expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        protein_list
            Return protein expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        library_size
            Scale the expression frequencies to a common library size.
            This allows gene expression levels to be interpreted on a common scale of relevant
            magnitude.
        n_samples
            Get sample scale from multiple samples.
        sample_protein_mixing
            Sample mixing bernoulli, setting background to zero
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a `np.ndarray` instead of a `pd.DataFrame`. Includes gene
            names as columns. If either n_samples=1 or return_mean=True, defaults to False.
            Otherwise, it defaults to True.

        Returns
        -------
        - **gene_normalized_expression** - normalized expression for RNA
        - **protein_normalized_expression** - normalized expression for proteins

        If ``n_samples`` > 1 and ``return_mean`` is False, then the shape is ``(samples, cells, genes)``.
        Otherwise, shape is ``(cells, genes)``. Return type is ``pd.DataFrame`` unless ``return_numpy`` is True.

        """

        adata = adata if adata is not None else self.adata
        post = self._make_posterior(adata=adata, indices=indices)

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]
        if protein_list is None:
            protein_mask = slice(None)
        else:
            all_proteins = adata.uns["scvi_protein_names"]
            protein_mask = [True if p in protein_list else False for p in all_proteins]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                logger.warning(
                    "return_numpy must be True if n_samples > 1 and return_mean is False, returning np.ndarray"
                )
            return_numpy = True

        scale_list_gene = []
        scale_list_pro = []
        if (transform_batch is None) or (isinstance(transform_batch, int)):
            transform_batch = [transform_batch]
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
            batch_index = tensors[_CONSTANTS.BATCH_KEY]
            label = tensors[_CONSTANTS.LABELS_KEY]
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
                if library_size == "observed":
                    px_scale += outputs["px_"]["rate"]
                else:
                    px_scale += outputs["px_"]["scale"]
                px_scale = px_scale[..., gene_mask]

                py_ = outputs["py_"]
                # probability of background
                protein_mixing = 1 / (1 + torch.exp(-py_["mixing"]))
                if sample_protein_mixing is True:
                    protein_mixing = torch.distributions.Bernoulli(
                        protein_mixing
                    ).sample()
                    protein_mixing = protein_mixing[..., protein_mask]
                py_scale += py_["rate_fore"][..., protein_mask] * (1 - protein_mixing)
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

        if return_mean is True and n_samples > 1:
            scale_list_gene = torch.mean(scale_list_gene, dim=-1)
            scale_list_pro = torch.mean(scale_list_pro, dim=-1)

        scale_list_gene = scale_list_gene.cpu().numpy()
        scale_list_pro = scale_list_pro.cpu().numpy()

        if return_numpy is None or return_numpy is False:
            gene_df = pd.DataFrame(
                scale_list_gene,
                columns=adata.var_names[gene_mask],
                index=adata.obs_names,
            )
            pro_df = pd.DataFrame(
                scale_list_pro,
                columns=adata.uns["scvi_protein_names"][protein_mask],
                index=adata.obs_names,
            )

            return gene_df, pro_df
        else:
            return scale_list_gene, scale_list_pro

    def differential_expression(
        self, groupby, group1=None, group2="rest", adata=None, within_key=None
    ):
        # group 1 and group 2 are valid obs keys in the anndata
        # runs 1vsall or 1v1 based on group1 and group2 choices
        # runs within cluster
        # new differential expression class

        raise NotImplementedError

        # if group2 == "rest":
        #     idx2 = ~group1
        # else:
        #     idx2 = 0

        # DE = DifferentialExpression(self.get_sample_scale)
        # pass

    @torch.no_grad()
    def posterior_predictive_sample(
        self,
        adata=None,
        indices=None,
        n_samples: int = 1,
        gene_list: Union[list, np.ndarray] = None,
        protein_list: Union[list, np.ndarray] = None,
    ) -> np.ndarray:
        r"""Generate observation samples from the posterior predictive distribution

        The posterior predictive distribution is written as :math:`p(\hat{x}, \hat{y} \mid x, y)`.

        Parameters
        ----------
        n_samples
            Number of required samples for each cell
        gene_list
            Indices or names of genes of interest

        Returns
        -------
        x_new : :py:class:`np.ndarray`
            tensor with shape (n_cells, n_genes, n_samples)
        """
        assert self.model.reconstruction_loss_gene in ["nb"]

        adata = adata if adata is not None else self.adata
        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]
        if protein_list is None:
            protein_mask = slice(None)
        else:
            all_proteins = adata.uns["scvi_protein_names"]
            protein_mask = [True if p in protein_list else False for p in all_proteins]

        post = self._make_posterior(adata=adata, indices=indices)

        posterior_list = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            batch_idx = tensors[_CONSTANTS.BATCH_KEY]
            labels = tensors[_CONSTANTS.LABELS_KEY]
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
            with torch.no_grad():
                outputs = self.model.inference(
                    x, y, batch_index=batch_idx, label=labels, n_samples=n_samples
                )
            px_ = outputs["px_"]
            py_ = outputs["py_"]

            pi = 1 / (1 + torch.exp(-py_["mixing"]))
            mixing_sample = torch.distributions.Bernoulli(pi).sample()
            protein_rate = (
                py_["rate_fore"] * (1 - mixing_sample)
                + py_["rate_back"] * mixing_sample
            )
            rate = torch.cat(
                (px_["rate"][..., gene_mask], protein_rate[..., protein_mask]), dim=-1
            )
            if len(px_["r"].size()) == 2:
                px_dispersion = px_["r"]
            else:
                px_dispersion = torch.ones_like(x) * px_["r"]
            if len(py_["r"].size()) == 2:
                py_dispersion = py_["r"]
            else:
                py_dispersion = torch.ones_like(y) * py_["r"]

            dispersion = torch.cat(
                (px_dispersion[..., gene_mask], py_dispersion[..., protein_mask]),
                dim=-1,
            )

            # This gamma is really l*w using scVI manuscript notation
            p = rate / (rate + dispersion)
            r = dispersion
            l_train = torch.distributions.Gamma(r, (1 - p) / p).sample()
            data = torch.distributions.Poisson(l_train).sample().cpu().numpy()
            # """
            # In numpy (shape, scale) => (concentration, rate), with scale = p /(1 - p)
            # rate = (1 - p) / p  # = 1/scale # used in pytorch
            # """
            posterior_list += [data]

            posterior_list[-1] = np.transpose(posterior_list[-1], (1, 2, 0))

        posterior_list = np.concatenate(posterior_list, axis=0)

        return posterior_list

    @torch.no_grad()
    def _get_denoised_samples(
        self,
        adata=None,
        indices=None,
        n_samples: int = 25,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[int] = None,
    ) -> np.ndarray:
        """Return samples from an adjusted posterior predictive.

        Parameters
        ----------
        n_samples
            How may samples per cell
        batch_size
            Mini-batch size for sampling. Lower means less GPU memory footprint
        rna_size_factor
            size factor for RNA prior to sampling gamma distribution
        transform_batch
            int of which batch to condition on for all cells

        Returns
        -------

        """

        post = self._make_posterior(adata=adata, indices=indices)

        posterior_list = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            batch_idx = tensors[_CONSTANTS.BATCH_KEY]
            label = tensors[_CONSTANTS.LABELS_KEY]
            y = tensors[_CONSTANTS.PROTEIN_EXP_KEY]
            with torch.no_grad():
                outputs = self.model.inference(
                    x,
                    y,
                    batch_index=batch_idx,
                    label=label,
                    n_samples=n_samples,
                    transform_batch=transform_batch,
                )
            px_ = outputs["px_"]
            py_ = outputs["py_"]

            pi = 1 / (1 + torch.exp(-py_["mixing"]))
            mixing_sample = torch.distributions.Bernoulli(pi).sample()
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
            l_train = torch.distributions.Gamma(r, (1 - p) / p).sample()
            data = l_train.cpu().numpy()
            # make background 0
            data[:, :, self.adata.shape[1] :] = (
                data[:, :, self.adata.shape[1] :] * (1 - mixing_sample).cpu().numpy()
            )
            posterior_list += [data]

            posterior_list[-1] = np.transpose(posterior_list[-1], (1, 2, 0))

        return np.concatenate(posterior_list, axis=0)

    @torch.no_grad()
    def get_feature_correlation_matrix(
        self,
        adata=None,
        indices=None,
        n_samples: int = 10,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[Union[int, List[int]]] = None,
        correlation_type: Literal["spearman", "pearson"] = "spearman",
        log_transform: bool = False,
    ) -> pd.DataFrame:
        """Generate gene-gene correlation matrix using scvi uncertainty and expression

        Parameters
        ----------
        n_samples
            How may samples per cell
        batch_size
            Mini-batch size for sampling. Lower means less GPU memory footprint
        rna_size_factor
            size factor for RNA prior to sampling gamma distribution
        transform_batch
            Batches to condition on.
            If transform_batch is:

            - None, then real observed batch is used
            - int, then batch transform_batch is used
            - list of int, then values are averaged over provided batches.
        correlation_type
            One of "pearson", "spearman"
        log_transform
            Whether to log transform denoised values prior to correlation calculation

        Returns
        -------
        Gene-protein-gene-protein correlation matrix
        """

        from scipy.stats import spearmanr

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
                flattened[:, : self.gene_dataset.n_genes] = np.log(
                    flattened[:, : self.gene_dataset.n_genes] + 1e-8
                )
                flattened[:, self.gene_dataset.n_genes :] = np.log1p(
                    flattened[:, self.gene_dataset.n_genes :]
                )
            if correlation_type == "pearson":
                corr_matrix = np.corrcoef(flattened, rowvar=False)
            else:
                corr_matrix = spearmanr(flattened, axis=0)
            corr_mats.append(corr_matrix)
        corr_matrix = np.mean(np.stack(corr_mats), axis=0)

        names = np.concatenate(
            [np.asarray(self.adata.var_names), self.adata.uns["scvi_protein_names"]]
        )
        return pd.DataFrame(corr_matrix, index=names, columns=names)

    # TODO
    @torch.no_grad()
    def get_likelihood_parameters(
        self,
        adata=None,
        indices=None,
        n_samples: Optional[int] = 1,
        give_mean: Optional[bool] = False,
    ) -> Dict[str, np.ndarray]:

        r"""Estimates for the parameters of the likelihood :math:`p(x \mid z)`.


        """

        post = self._make_posterior(adata=adata, indices=indices)

        rna_dropout_list = []
        rna_mean_list = []
        rna_dispersion_list = []
        for tensors in post:
            x = tensors[_CONSTANTS.X_KEY]
            batch_idx = tensors[_CONSTANTS.BATCH_KEY]
            labels = tensors[_CONSTANTS.LABELS_KEY]

            outputs = self.model.inference(
                x, batch_index=batch_idx, y=labels, n_samples=n_samples
            )
            px_r = outputs["px_r"]
            px_rate = outputs["px_rate"]
            px_dropout = outputs["px_dropout"]

            n_batch = px_rate.size(0) if n_samples == 1 else px_rate.size(1)
            rna_dispersion_list += [
                np.repeat(np.array(px_r.cpu())[np.newaxis, :], n_batch, axis=0)
            ]
            rna_mean_list += [np.array(px_rate.cpu())]
            rna_dropout_list += [np.array(px_dropout.cpu())]

        dropout = np.concatenate(rna_dropout_list)
        means = np.concatenate(rna_mean_list)
        dispersions = np.concatenate(rna_dispersion_list)
        if give_mean and n_samples > 1:
            dropout = dropout.mean(0)
            means = means.mean(0)

        return_dict = {}
        return_dict["RNA_mean"] = means

        if self.model.reconstruction_loss_gene == "zinb":
            return_dict["RNA_dropout"] = dropout
            return_dict["RNA_dispersions"] = dispersions
        if self.model.reconstruction_loss_gene == "nb":
            return_dict["RNA_dispersions"] = dispersions

        return return_dict
