import logging
from typing import Optional

import numpy as np
import pandas as pd
import pyro
import torch
from anndata import AnnData

from scvi._constants import _CONSTANTS
from scvi.module import LDAPyroModule

from .base import BaseModelClass, PyroSviTrainMixin

logger = logging.getLogger(__name__)


class LDA(PyroSviTrainMixin, BaseModelClass):
    """
    Latent Dirichlet Allocation [Blei03]_.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    n_components
        Number of components/topics to model.
    n_hidden
        Number of nodes in the hidden layer of the encoder.
    cell_component_prior
        Prior of cell component distribution. If `None`, defaults to `1 / n_components`.
    component_gene_prior
        Prior of component gene distribution. If `None`, defaults to `1 / n_components`.
    """

    def __init__(
        self,
        adata: AnnData,
        n_components: int = 20,
        n_hidden: int = 128,
        cell_component_prior: Optional[float] = None,
        component_gene_prior: Optional[float] = None,
    ):
        # in case any other model was created before that shares the same parameter names.
        pyro.clear_param_store()

        super().__init__(adata)

        self.module = LDAPyroModule(
            n_input=self.summary_stats["n_vars"],
            n_components=n_components,
            n_hidden=n_hidden,
            cell_component_prior=cell_component_prior,
            component_gene_prior=component_gene_prior,
        )

    def _check_var_equality(self, adata: AnnData):
        """
        Checks that the var names are equivalent between the source AnnData object and
        the AnnData passed in. Throws a ValueError if not.

        Parameters
        ----------
        adata
            AnnData to compare against.
        """
        source_var_names = self.adata.var_names.astype(str)
        user_var_names = adata.var_names.astype(str)
        if not np.array_equal(user_var_names, source_var_names):
            raise ValueError(
                "`adata` passed into `transform` does not have matching var_names "
                "with the source adata the model was trained with."
            )

    def _check_if_not_trained(self):
        """
        Check if the model is not trained. Throws an AssertionError if not.
        """
        assert (
            self.is_trained_
        ), "Trying to query inferred values from an untrained model. Please train the model first."

    def get_components(self) -> pd.DataFrame:
        """
        Gets the component to gene transition matrix.

        Parameters
        ----------
        adata
            AnnData to transform. If None, returns the component to gene transition matrix for
            the source AnnData.

        Returns
        -------
        A `n_components x n_var` Pandas DataFrame containing the component to gene transition matrix.
        """
        self._check_if_not_trained()

        return pd.DataFrame(
            data=self.module.components.numpy(), columns=self.adata.var_names
        )

    def transform(self, adata: Optional[AnnData] = None) -> pd.DataFrame:
        """
        Converts a count matrix to an inferred component distribution.

        Parameters
        ----------
        adata
            AnnData to transform. If None, returns the component distribution for the source AnnData.

        Returns
        -------
        A `n_obs x n_components` Pandas DataFrame containing the normalized estimate
        of the component distribution for each observation.
        """
        if adata is not None:
            self._check_var_equality(adata)
        self._check_if_not_trained()

        user_adata = adata or self.adata
        dl = self._make_data_loader(
            adata=user_adata, indices=np.arange(user_adata.n_obs)
        )

        transformed_xs = []
        for tensors in dl:
            x = tensors[_CONSTANTS.X_KEY]
            transformed_xs.append(self.module.transform(x))

        transformed_x = torch.cat(transformed_xs).numpy()
        return pd.DataFrame(data=transformed_x, index=user_adata.obs_names)

    def perplexity(self, adata: Optional[AnnData] = None) -> float:
        """
        Computes the approximate perplexity of the for `adata`, where perplexity is defined
        as exp(-1 * log-likelihood per count).

        Parameters
        ----------
        adata
            AnnData to compute perplexity for. If None, returns the perplexity for the source AnnData.

        Returns
        -------
        Perplexity as a float.
        """
        if adata is not None:
            self._check_var_equality(adata)
        self._check_if_not_trained()

        user_adata = adata or self.adata
        dl = self._make_data_loader(
            adata=user_adata, indices=np.arange(user_adata.n_obs)
        )

        perplexities = []
        batch_counts = []
        for tensors in dl:
            x = tensors[_CONSTANTS.X_KEY]
            x_counts = x.sum().item()
            batch_counts.append(x_counts)
            perplexities.append(self.module.perplexity(x))

        normalized_batch_counts = np.array(batch_counts) / np.sum(batch_counts)
        return np.prod(np.power(perplexities, normalized_batch_counts))
