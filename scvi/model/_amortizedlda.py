import collections.abc
import logging
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
import pyro
import torch
from anndata import AnnData

from scvi._constants import _CONSTANTS
from scvi._docs import dsp
from scvi.data._anndata import _setup_anndata
from scvi.module import AmortizedLDAPyroModule

from .base import BaseModelClass, PyroSviTrainMixin

logger = logging.getLogger(__name__)


class AmortizedLDA(PyroSviTrainMixin, BaseModelClass):
    """
    Amortized Latent Dirichlet Allocation [Blei03]_.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.AmortizedLDA.setup_anndata`.
    n_topics
        Number of topics to model.
    n_hidden
        Number of nodes in the hidden layer of the encoder.
    cell_topic_prior
        Prior of cell topic distribution. If `None`, defaults to `1 / n_topics`.
    topic_gene_prior
        Prior of topic gene distribution. If `None`, defaults to `1 / n_topics`.

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.AmortizedLDA.setup_anndata(adata)
    >>> model = scvi.model.AmortizedLDA(adata)
    >>> model.train()
    >>> gene_by_topic = model.get_gene_by_topic()
    >>> adata.obsm["X_LDA"] = model.get_latent_representation()
    """

    def __init__(
        self,
        adata: AnnData,
        n_topics: int = 20,
        n_hidden: int = 128,
        cell_topic_prior: Optional[Union[float, Sequence[float]]] = None,
        topic_gene_prior: Optional[Union[float, Sequence[float]]] = None,
    ):
        # in case any other model was created before that shares the same parameter names.
        pyro.clear_param_store()

        super().__init__(adata)

        n_input = self.summary_stats["n_vars"]

        if (
            cell_topic_prior is not None
            and not isinstance(cell_topic_prior, float)
            and (
                not isinstance(cell_topic_prior, collections.abc.Sequence)
                or len(cell_topic_prior) != n_topics
            )
        ):
            raise ValueError(
                f"cell_topic_prior, {cell_topic_prior}, must be None, "
                f"a float or a Sequence of length n_topics."
            )
        if (
            topic_gene_prior is not None
            and not isinstance(topic_gene_prior, float)
            and (
                not isinstance(topic_gene_prior, collections.abc.Sequence)
                or len(topic_gene_prior) != n_input
            )
        ):
            raise ValueError(
                f"topic_gene_prior, {topic_gene_prior}, must be None, "
                f"a float or a Sequence of length n_input."
            )

        self.module = AmortizedLDAPyroModule(
            n_input=n_input,
            n_topics=n_topics,
            n_hidden=n_hidden,
            cell_topic_prior=cell_topic_prior,
            topic_gene_prior=topic_gene_prior,
        )

        self.init_params_ = self._get_init_params(locals())

    @staticmethod
    @dsp.dedent
    def setup_anndata(
        adata: AnnData,
        layer: Optional[str] = None,
        copy: bool = False,
    ) -> Optional[AnnData]:
        """
        %(setup_anndata_summary)s

        Parameters
        ----------
        %(setup_anndata_param_adata)s
        %(setup_anndata_param_layer)s
        %(setup_anndata_param_copy)s

        Returns
        -------
        %(setup_anndata_returns)s
        """
        return _setup_anndata(
            adata,
            layer=layer,
            copy=copy,
        )

    def get_gene_by_topic(self, give_mean=True) -> pd.DataFrame:
        """
        Gets the gene by topic matrix.

        Parameters
        ----------
        adata
            AnnData to transform. If None, returns the gene by topic matrix for
            the source AnnData.
        give_mean
            Give mean of distribution if True or sample from it.

        Returns
        -------
        A `n_var x n_topics` Pandas DataFrame containing the gene by topic matrix.
        """
        self._check_if_trained(warn=False)

        topic_by_gene = self.module.topic_by_gene(give_mean=give_mean)

        return pd.DataFrame(
            data=topic_by_gene.numpy().T,
            index=self.adata.var_names,
            columns=[f"topic_{i}" for i in range(topic_by_gene.shape[0])],
        )

    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        batch_size: Optional[int] = None,
        give_mean: bool = True,
    ) -> pd.DataFrame:
        """
        Converts a count matrix to an inferred topic distribution.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        give_mean
            Give mean of distribution or sample from it.

        Returns
        -------
        A `n_obs x n_topics` Pandas DataFrame containing the normalized estimate
        of the topic distribution for each observation.
        """
        self._check_if_trained(warn=False)
        adata = self._validate_anndata(adata)

        dl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        transformed_xs = []
        for tensors in dl:
            x = tensors[_CONSTANTS.X_KEY]
            transformed_xs.append(
                self.module.get_topic_distribution(x, give_mean=give_mean)
            )
        transformed_x = torch.cat(transformed_xs).numpy()

        return pd.DataFrame(
            data=transformed_x,
            index=adata.obs_names,
            columns=[f"topic_{i}" for i in range(transformed_x.shape[1])],
        )

    def get_elbo(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        batch_size: Optional[int] = None,
    ) -> float:
        """
        Return the ELBO for the data.

        The ELBO is a lower bound on the log likelihood of the data used for optimization
        of VAEs. Note, this is not the negative ELBO, higher is better.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        The positive ELBO.
        """
        self._check_if_trained(warn=False)
        adata = self._validate_anndata(adata)

        dl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        elbos = []
        for tensors in dl:
            x = tensors[_CONSTANTS.X_KEY]
            library = x.sum(dim=1)
            elbos.append(self.module.get_elbo(x, library, len(dl.indices)))
        return np.mean(elbos)

    def get_perplexity(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        batch_size: Optional[int] = None,
    ) -> float:
        """
        Computes approximate perplexity for `adata`.

        Perplexity is defined as exp(-1 * log-likelihood per count).

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        Perplexity.
        """
        self._check_if_trained(warn=False)
        adata = self._validate_anndata(adata)

        dl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        total_counts = sum(tensors[_CONSTANTS.X_KEY].sum().item() for tensors in dl)

        return np.exp(
            self.get_elbo(adata=adata, indices=indices, batch_size=batch_size)
            / total_counts
        )
