import copy
import logging
import warnings
from typing import List, Optional, Union, Tuple, Dict

import numpy as np
import torch
import anndata
import torch.distributions as distributions
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score as ARI
from sklearn.metrics import normalized_mutual_info_score as NMI
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture as GMM
from torch.utils.data import DataLoader

from scvi.dataset._biodataset import BioDataset
from scvi.inference.posterior_utils import (
    entropy_batch_mixing,
    nn_overlap,
    unsupervised_clustering_accuracy,
    knn_purity,
)
from scvi.models.log_likelihood import (
    compute_elbo,
    compute_reconstruction_error,
    compute_marginal_log_likelihood_scvi,
    compute_marginal_log_likelihood_autozi,
)
from scvi.models.distributions import NegativeBinomial, ZeroInflatedNegativeBinomial
from scipy.stats import spearmanr
from scvi import _CONSTANTS

logger = logging.getLogger(__name__)


class BatchSampler(torch.utils.data.sampler.Sampler):
    """Custom torch Sampler that returns a list of indices of size batch_size

    Parameters
    ----------
    indices
        list of indices to sample from
    batch_size
        batch size of each iteration
    shuffle
        if ``True``, shuffles indices before sampling
    """

    def __init__(self, indices: np.ndarray, batch_size: int, shuffle: bool):
        self.indices = indices
        self.batch_size = batch_size
        self.shuffle = shuffle

    def __iter__(self):
        if self.shuffle is True:
            idx = torch.randperm(len(self.indices)).tolist()
        else:
            idx = torch.arange(len(self.indices)).tolist()

        data_iter = iter(
            [
                self.indices[idx[i : i + self.batch_size]]
                for i in range(0, len(idx), self.batch_size)
            ]
        )
        return data_iter

    def __len__(self):
        return len(self.indices) // self.batch_size


class Posterior:
    """The functional data unit.

    A `Posterior` instance is instantiated with a model and a gene_dataset, and
    as well as additional arguments that for Pytorch's `DataLoader`. A subset of indices can be specified, for
    purposes such as splitting the data into train/test or labelled/unlabelled (for semi-supervised learning).
    Each trainer instance of the `Trainer` class can therefore have multiple `Posterior` instances to train a model.
    A `Posterior` instance also comes with many methods or utilities for its corresponding data.

    Parameters
    ----------
    model
        A model instance from class ``VAE``, ``VAEC``, ``SCANVI``
    gene_dataset
        A gene_dataset instance like ``CortexDataset()``
    shuffle
        Specifies if a `RandomSampler` or a `SequentialSampler` should be used
    indices
        Specifies how the data should be split with regards to train/test or labelled/unlabelled
    use_cuda
        Default: ``True``
    data_loader_kwargs
        Keyword arguments to passed into the `DataLoader`

    Examples
    --------
    Let us instantiate a `trainer`, with a gene_dataset and a model

    A `UnsupervisedTrainer` instance has two `Posterior` attributes: `train_set` and `test_set`
    For this subset of the original gene_dataset instance, we can examine the differential expression,
    log_likelihood, entropy batch mixing, etc.

    >>> gene_dataset = CortexDataset()
    >>> vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
    ... n_labels=gene_dataset.n_labels, use_cuda=True)
    >>> trainer = UnsupervisedTrainer(vae, gene_dataset)
    >>> trainer.train(n_epochs=50)

    >>> trainer.train_set.differential_expression_stats()
    >>> trainer.train_set.reconstruction_error()
    >>> trainer.train_set.entropy_batch_mixing()
    >>> trainer.train_set.show_t_sne(n_samples=1000, color_by="labels")

    """

    def __init__(
        self,
        model,
        adata: anndata.AnnData,
        shuffle=False,
        indices=None,
        use_cuda=True,
        batch_size=128,
        data_loader_kwargs=dict(),
    ):
        self.model = model
        assert "scvi_data_registry" in adata.uns.keys(), ValueError(
            "Please run setup_anndata() on your anndata object first."
        )
        for key in self._data_and_attributes.keys():
            assert key in adata.uns["scvi_data_registry"].keys(), ValueError(
                "{} required for model but not included when setup_anndata was run".format(
                    key
                )
            )
        self.gene_dataset = BioDataset(adata, getitem_tensors=self._data_and_attributes)
        self.to_monitor = []
        self.use_cuda = use_cuda

        if indices is None:
            inds = np.arange(len(self.gene_dataset))
            if shuffle:
                sampler_kwargs = {
                    "indices": inds,
                    "batch_size": batch_size,
                    "shuffle": True,
                }
            else:
                sampler_kwargs = {
                    "indices": inds,
                    "batch_size": batch_size,
                    "shuffle": False,
                }
        else:
            if hasattr(indices, "dtype") and indices.dtype is np.dtype("bool"):
                indices = np.where(indices)[0].ravel()
            indices = np.asarray(indices)
            sampler_kwargs = {
                "indices": indices,
                "batch_size": batch_size,
                "shuffle": True,
            }

        self.sampler_kwargs = sampler_kwargs
        sampler = BatchSampler(**self.sampler_kwargs)
        self.data_loader_kwargs = copy.copy(data_loader_kwargs)
        # do not touch batch size here, sampler gives batched indices
        self.data_loader_kwargs.update({"sampler": sampler, "batch_size": None})
        self.data_loader = DataLoader(self.gene_dataset, **self.data_loader_kwargs)
        self.original_indices = self.indices

    @property
    def _data_and_attributes(self):
        return {
            _CONSTANTS.X_KEY: np.float32,
            _CONSTANTS.BATCH_KEY: np.int64,
            _CONSTANTS.LOCAL_L_MEAN_KEY: np.float32,
            _CONSTANTS.LOCAL_L_VAR_KEY: np.float32,
            _CONSTANTS.LABELS_KEY: np.int64,
        }

    def accuracy(self):
        pass

    accuracy.mode = "max"

    # def save_posterior(self, dir_path: str):
    #     """Saves the posterior properties in folder `dir_path`.

    #     To ensure safety, this method requires that `dir_path` does not exist.
    #     The posterior can then be retrieved later on with the function `load_posterior`

    #     Parameters
    #     ----------
    #     dir_path
    #         non-existing directory in which the posterior properties will be saved.
    #     """

    #     if not os.path.exists(dir_path):
    #         os.makedirs(dir_path)
    #     else:
    #         raise ValueError(
    #             "{} already exists. Please provide an unexisting directory for saving.".format(
    #                 dir_path
    #             )
    #         )
    #     anndata_dataset = self.gene_dataset.to_anndata()

    #     anndata_dataset.write(
    #         os.path.join(dir_path, "anndata_dataset.h5ad"), compression="lzf"
    #     )
    #     with open(os.path.join(dir_path, "posterior_type.txt"), "w") as post_file:
    #         post_file.write(self.posterior_type)
    #     torch.save(self.model.state_dict(), os.path.join(dir_path, "model_params.pt"))

    #     # Saves posterior indices and kwargs that can easily be retrieved
    #     data_loader_kwargs = pd.Series(
    #         {key: vals for key, vals in self.data_loader_kwargs.items()}
    #     )
    #     data_loader_kwargs = data_loader_kwargs[
    #         ~data_loader_kwargs.index.isin(["collate_fn", "sampler"])
    #     ]
    #     data_loader_kwargs.to_hdf(
    #         os.path.join(dir_path, "data_loader_kwargs.h5"), key="data_loader"
    #     )
    #     np.save(file=os.path.join(dir_path, "indices.npy"), arr=np.array(self.indices))

    #     sampler_kwargs = pd.Series(
    #         {key: vals for key, vals in self.sampler_kwargs.items()}
    #     )
    #     sampler_kwargs = sampler_kwargs[~sampler_kwargs.index.isin(["indices"])]
    #     sampler_kwargs.to_hdf(
    #         os.path.join(dir_path, "sampler_kwargs.h5"), key="sampler"
    #     )

    @property
    def indices(self) -> np.ndarray:
        """Returns the current dataloader indices used by the object"""
        if hasattr(self.data_loader.sampler, "indices"):
            return self.data_loader.sampler.indices
        else:
            return np.arange(len(self.gene_dataset))

    @property
    def n_cells(self) -> int:
        """returns the number of studied cells."""
        if hasattr(self.data_loader.sampler, "indices"):
            return len(self.data_loader.sampler.indices)
        else:
            return self.gene_dataset.n_cells

    @property
    def posterior_type(self) -> str:
        """Returns the posterior class name"""
        return self.__class__.__name__

    def __iter__(self):
        return map(self.to_cuda, iter(self.data_loader))

    def to_cuda(self, tensors: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
        """Converts list of tensors to cuda.

        Parameters
        ----------
        tensors
            tensors to convert
        """
        return {k: (t.cuda() if self.use_cuda else t) for k, t in tensors.items()}

    def update(self, data_loader_kwargs: dict) -> "Posterior":
        """Updates the dataloader

        Parameters
        ----------
        data_loader_kwargs
            dataloader updates.

        Returns
        -------
        Updated posterior
        """
        posterior = copy.copy(self)
        posterior.data_loader_kwargs = copy.copy(self.data_loader_kwargs)
        posterior.data_loader_kwargs.update(data_loader_kwargs)
        posterior.data_loader = DataLoader(
            self.gene_dataset, **posterior.data_loader_kwargs
        )
        return posterior

    def update_batch_size(self, batch_size):
        self.sampler_kwargs.update({"batch_size": batch_size})
        sampler = BatchSampler(**self.sampler_kwargs)
        return self.update({"sampler": sampler, "batch_size": None})

    def sequential(self, batch_size: Optional[int] = 128) -> "Posterior":
        """Returns a copy of the object that iterate over the data sequentially.

        Parameters
        ----------
        batch_size
            New batch size.

        """
        self.sampler_kwargs = {
            "indices": self.indices,
            "batch_size": batch_size,
            "shuffle": False,
        }
        return self.update({"sampler": BatchSampler(**self.sampler_kwargs)})

    @torch.no_grad()
    def elbo(self) -> torch.Tensor:
        """Returns the Evidence Lower Bound associated to the object."""
        elbo = compute_elbo(self.model, self)
        logger.debug("ELBO : %.4f" % elbo)
        return elbo

    elbo.mode = "min"

    @torch.no_grad()
    def reconstruction_error(self) -> torch.Tensor:
        """Returns the reconstruction error associated to the object."""
        reconstruction_error = compute_reconstruction_error(self.model, self)
        logger.debug("Reconstruction Error : %.4f" % reconstruction_error)
        return reconstruction_error

    reconstruction_error.mode = "min"

    @torch.no_grad()
    def marginal_ll(self, n_mc_samples: Optional[int] = 1000) -> torch.Tensor:
        """Estimates the marginal likelihood of the object's data.

        Parameters
        ----------
        n_mc_samples
            Number of MC estimates to use

        Returns
        -------
        Marginal LL
        """
        if (
            hasattr(self.model, "reconstruction_loss")
            and self.model.reconstruction_loss == "autozinb"
        ):
            ll = compute_marginal_log_likelihood_autozi(self.model, self, n_mc_samples)
        else:
            ll = compute_marginal_log_likelihood_scvi(self.model, self, n_mc_samples)
        logger.debug("True LL : %.4f" % ll)
        return ll

    @torch.no_grad()
    def entropy_batch_mixing(self, **kwargs) -> torch.Tensor:
        """Returns the object's entropy batch mixing.
        """
        if self.gene_dataset.n_batches == 2:
            latent, batch_indices, labels = self.get_latent()
            be_score = entropy_batch_mixing(latent, batch_indices, **kwargs)
            logger.debug("Entropy batch mixing : {}".format(be_score))
            return be_score

    entropy_batch_mixing.mode = "max"

    def update_sampler_indices(self, idx: Union[List, np.ndarray]):
        """Updates the dataloader indices.

        More precisely, this method can be used to temporarily change which cells __iter__
        will yield. This is particularly useful for computational considerations when one is only interested
        in a subset of the cells of the Posterior object.
        This method should be used carefully and requires to reset the dataloader to its
        original value after use.

        Parameters
        ----------
        idx :
            Indices (in [0, len(dataset)] to sample from

        Examples
        --------
        >>> old_loader = self.data_loader
        >>> cell_indices = np.array([1, 2, 3])
        >>> self.update_sampler_indices(cell_indices)
        >>> for tensors in self:
        >>>    # your code

        >>> # Do not forget next line!
        >>> self.data_loader = old_loader
        """
        self.sampler_kwargs.update({"indices": idx})
        sampler = BatchSampler(**self.sampler_kwargs)
        self.data_loader_kwargs.update({"sampler": sampler, "batch_size": None})
        self.data_loader = DataLoader(self.gene_dataset, **self.data_loader_kwargs)

    @torch.no_grad()
    def differential_expression_stats(self, M_sampling: int = 100) -> Tuple:
        """Output average over statistics in a symmetric way (a against b), forget the sets if permutation is True

        Parameters
        ----------
        M_sampling
            number of samples

        Returns
        -------
        type
            Tuple px_scales, all_labels where (i) px_scales: scales of shape (M_sampling, n_genes)
            (ii) all_labels: labels of shape (M_sampling, )

        """

        warnings.warn(
            "differential_expression_stats() is deprecated; "
            "use differential_expression_score() or get_sample_scale().",
            category=DeprecationWarning,
        )

        px_scales = []
        all_labels = []

        batch_size = max(
            self.sampler_kwargs["batch_size"] // M_sampling, 2
        )  # Reduce batch_size on GPU
        if len(self.gene_dataset) % batch_size == 1:
            batch_size += 1
        for tensors in self.update_batch_size(batch_size):
            sample_batch, _, _, batch_index, labels = self._unpack_tensors(tensors)
            px_scales += [
                np.array(
                    (
                        self.model.get_sample_scale(
                            sample_batch,
                            batch_index=batch_index,
                            y=labels,
                            n_samples=M_sampling,
                        )
                    ).cpu()
                )
            ]

            # Align the sampling
            if M_sampling > 1:
                px_scales[-1] = (px_scales[-1].transpose((1, 0, 2))).reshape(
                    -1, px_scales[-1].shape[-1]
                )
            all_labels += [np.array((labels.repeat(1, M_sampling).view(-1, 1)).cpu())]

        px_scales = np.concatenate(px_scales)
        all_labels = np.concatenate(all_labels).ravel()  # this will be used as boolean

        return px_scales, all_labels

    @torch.no_grad()
    def generate(
        self,
        n_samples: int = 100,
        genes: Union[list, np.ndarray] = None,
        batch_size: int = 128,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Create observation samples from the Posterior Predictive distribution

        Parameters
        ----------
        n_samples
            Number of required samples for each cell
        genes
            Indices of genes of interest
        batch_size
            Desired Batch size to generate data

        Returns
        -------
        x_new : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_genes, n_samples)
        x_old : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_genes)

        """
        assert self.model.reconstruction_loss in ["zinb", "nb", "poisson"]
        x_old = []
        x_new = []
        for tensors in self.update_batch_size(batch_size):
            sample_batch, _, _, batch_index, labels = self._unpack_tensors(tensors)

            outputs = self.model.inference(
                sample_batch, batch_index=batch_index, y=labels, n_samples=n_samples
            )
            px_r = outputs["px_r"]
            px_rate = outputs["px_rate"]
            px_dropout = outputs["px_dropout"]

            if self.model.reconstruction_loss == "poisson":
                l_train = px_rate
                l_train = torch.clamp(l_train, max=1e8)
                dist = distributions.Poisson(
                    l_train
                )  # Shape : (n_samples, n_cells_batch, n_genes)
            elif self.model.reconstruction_loss == "nb":
                dist = NegativeBinomial(mu=px_rate, theta=px_r)
            elif self.model.reconstruction_loss == "zinb":
                dist = ZeroInflatedNegativeBinomial(
                    mu=px_rate, theta=px_r, zi_logits=px_dropout
                )
            else:
                raise ValueError(
                    "{} reconstruction error not handled right now".format(
                        self.model.reconstruction_loss
                    )
                )
            gene_expressions = dist.sample().permute(
                [1, 2, 0]
            )  # Shape : (n_cells_batch, n_genes, n_samples)

            x_old.append(sample_batch.cpu())
            x_new.append(gene_expressions.cpu())

        x_old = torch.cat(x_old)  # Shape (n_cells, n_genes)
        x_new = torch.cat(x_new)  # Shape (n_cells, n_genes, n_samples)
        if genes is not None:
            gene_ids = self.gene_dataset.genes_to_index(genes)
            x_new = x_new[:, gene_ids, :]
            x_old = x_old[:, gene_ids]
        return x_new.numpy(), x_old.numpy()

    def _unpack_tensors(self, tensors):
        x = tensors[_CONSTANTS.X_KEY]
        local_l_mean = tensors[_CONSTANTS.LOCAL_L_MEAN_KEY]
        local_l_var = tensors[_CONSTANTS.LOCAL_L_VAR_KEY]
        batch_index = tensors[_CONSTANTS.BATCH_KEY]
        y = tensors[_CONSTANTS.LABELS_KEY]
        return x, local_l_mean, local_l_var, batch_index, y

    @torch.no_grad()
    def generate_denoised_samples(
        self,
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
        posterior_list = []
        for tensors in self.update_batch_size(batch_size):
            sample_batch, _, _, batch_index, labels = self._unpack_tensors(tensors)

            outputs = self.model.inference(
                sample_batch, batch_index=batch_index, y=labels, n_samples=n_samples
            )
            px_scale = outputs["px_scale"]
            px_r = outputs["px_r"]

            rate = rna_size_factor * px_scale
            if len(px_r.size()) == 2:
                px_dispersion = px_r
            else:
                px_dispersion = torch.ones_like(sample_batch) * px_r

            # This gamma is really l*w using scVI manuscript notation
            p = rate / (rate + px_dispersion)
            r = px_dispersion
            l_train = distributions.Gamma(r, (1 - p) / p).sample()
            data = l_train.cpu().numpy()
            # """
            # In numpy (shape, scale) => (concentration, rate), with scale = p /(1 - p)
            # rate = (1 - p) / p  # = 1/scale # used in pytorch
            # """
            posterior_list += [data]

            posterior_list[-1] = np.transpose(posterior_list[-1], (1, 2, 0))

        return np.concatenate(posterior_list, axis=0)

    @torch.no_grad()
    def generate_feature_correlation_matrix(
        self,
        n_samples: int = 10,
        batch_size: int = 64,
        rna_size_factor: int = 1000,
        transform_batch: Optional[Union[int, List[int]]] = None,
        correlation_type: str = "spearman",
    ) -> np.ndarray:
        """Wrapper of `generate_denoised_samples()` to create a gene-gene corr matrix

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

        Returns
        -------
        Gene-gene correlation matrix
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
            if correlation_type == "pearson":
                corr_matrix = np.corrcoef(flattened, rowvar=False)
            elif correlation_type == "spearman":
                corr_matrix, _ = spearmanr(flattened)
            else:
                raise ValueError(
                    "Unknown correlation type. Choose one of 'spearman', 'pearson'."
                )
            corr_mats.append(corr_matrix)
        corr_matrix = np.mean(np.stack(corr_mats), axis=0)
        return corr_matrix

    @torch.no_grad()
    def generate_parameters(
        self, n_samples: Optional[int] = 1, give_mean: Optional[bool] = False
    ) -> Tuple:

        """Estimates data's count means, dispersions and dropout logits.
        """
        dropout_list = []
        mean_list = []
        dispersion_list = []
        for tensors in self.sequential(1000):
            sample_batch, _, _, batch_index, labels = self._unpack_tensors(tensors)

            outputs = self.model.inference(
                sample_batch, batch_index=batch_index, y=labels, n_samples=n_samples
            )
            px_r = outputs["px_r"]
            px_rate = outputs["px_rate"]
            px_dropout = outputs["px_dropout"]

            n_batch = px_rate.size(0) if n_samples == 1 else px_rate.size(1)
            dispersion_list += [
                np.repeat(np.array(px_r.cpu())[np.newaxis, :], n_batch, axis=0)
            ]
            mean_list += [np.array(px_rate.cpu())]
            dropout_list += [np.array(px_dropout.cpu())]

        dropout = np.concatenate(dropout_list)
        means = np.concatenate(mean_list)
        dispersions = np.concatenate(dispersion_list)
        if give_mean and n_samples > 1:
            dropout = dropout.mean(0)
            means = means.mean(0)

        return (dropout, means, dispersions)

    @torch.no_grad()
    def get_stats(self) -> np.ndarray:
        libraries = []
        for tensors in self.sequential(batch_size=128):
            x, _, _, batch_index, y = self._unpack_tensors(tensors)

            library = self.model.inference(x, batch_index, y)["library"]
            libraries += [np.array(library.cpu())]
        libraries = np.concatenate(libraries)
        return libraries.ravel()

    @torch.no_grad()
    def knn_purity(self) -> torch.Tensor:
        """Computes kNN purity as described in [Lopez18]_"""
        latent, _, labels = self.get_latent()
        score = knn_purity(latent, labels)
        logger.debug("KNN purity score : {}".format(score))
        return score

    knn_purity.mode = "max"

    @torch.no_grad()
    def clustering_scores(self, prediction_algorithm: str = "knn") -> Tuple:
        if self.gene_dataset.adata.uns["scvi_summary_stats"]["n_labels"] > 1:
            latent, _, labels = self.get_latent()
            if prediction_algorithm == "knn":
                labels_pred = KMeans(
                    self.gene_dataset.adata.uns["scvi_summary_stats"]["n_labels"],
                    n_init=200,
                ).fit_predict(latent)
            elif prediction_algorithm == "gmm":
                gmm = GMM(self.gene_dataset.adata.uns["scvi_summary_stats"]["n_labels"])
                gmm.fit(latent)
                labels_pred = gmm.predict(latent)

            asw_score = silhouette_score(latent, labels)
            nmi_score = NMI(labels, labels_pred)
            ari_score = ARI(labels, labels_pred)
            uca_score = unsupervised_clustering_accuracy(labels, labels_pred)[0]
            logger.debug(
                "Clustering Scores:\nSilhouette: %.4f\nNMI: %.4f\nARI: %.4f\nUCA: %.4f"
                % (asw_score, nmi_score, ari_score, uca_score)
            )
            return asw_score, nmi_score, ari_score, uca_score

    @torch.no_grad()
    def nn_overlap_score(self, **kwargs) -> Tuple:
        """Quantify how much the similarity between cells in the mRNA latent space resembles their similarity at the
        protein level.

        Compute the overlap fold enrichment between the protein and mRNA-based cell 100-nearest neighbor
        graph and the Spearman correlation of the adjacency matrices.

        Parameters
        ----------
        **kwargs


        Returns
        -------

        """
        if hasattr(self.gene_dataset, "protein_expression_clr"):
            latent, _, _ = self.sequential().get_latent()
            protein_data = self.gene_dataset.protein_expression_clr[self.indices]
            spearman_correlation, fold_enrichment = nn_overlap(
                latent, protein_data, **kwargs
            )
            logger.debug(
                "Overlap Scores:\nSpearman Correlation: %.4f\nFold Enrichment: %.4f"
                % (spearman_correlation, fold_enrichment)
            )
            return spearman_correlation, fold_enrichment

    def raw_data(self) -> Tuple:
        """Returns raw data for classification"""
        data = self.gene_dataset[self.indices]
        return (data[_CONSTANTS.X_KEY], data[_CONSTANTS.LABELS_KEY].ravel())
