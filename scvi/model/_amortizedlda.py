import collections.abc
import logging
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
import pyro
import torch
from anndata import AnnData

from scvi._constants import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import LayerField
from scvi.module import AmortizedLDAPyroModule
from scvi.utils import setup_anndata_dsp

from .base import BaseModelClass, PyroSviTrainMixin

logger = logging.getLogger(__name__)


class AmortizedLDA(PyroSviTrainMixin, BaseModelClass):
    """
    Amortized Latent Dirichlet Allocation :cite:p:`Blei03`.

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
    topic_feature_prior
        Prior of topic feature distribution. If `None`, defaults to `1 / n_topics`.

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.AmortizedLDA.setup_anndata(adata)
    >>> model = scvi.model.AmortizedLDA(adata)
    >>> model.train()
    >>> feature_by_topic = model.get_feature_by_topic()
    >>> adata.obsm["X_LDA"] = model.get_latent_representation()
    """

    _module_cls = AmortizedLDAPyroModule

    def __init__(
        self,
        adata: AnnData,
        n_topics: int = 20,
        n_hidden: int = 128,
        cell_topic_prior: Optional[Union[float, Sequence[float]]] = None,
        topic_feature_prior: Optional[Union[float, Sequence[float]]] = None,
    ):
        # in case any other model was created before that shares the same parameter names.
        pyro.clear_param_store()

        super().__init__(adata)

        n_input = self.summary_stats.n_vars

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
            topic_feature_prior is not None
            and not isinstance(topic_feature_prior, float)
            and (
                not isinstance(topic_feature_prior, collections.abc.Sequence)
                or len(topic_feature_prior) != n_input
            )
        ):
            raise ValueError(
                f"topic_feature_prior, {topic_feature_prior}, must be None, "
                f"a float or a Sequence of length n_input."
            )

        self.module = self._module_cls(
            n_input=n_input,
            n_topics=n_topics,
            n_hidden=n_hidden,
            cell_topic_prior=cell_topic_prior,
            topic_feature_prior=topic_feature_prior,
        )

        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: Optional[str] = None,
        **kwargs,
    ) -> Optional[AnnData]:
        """
        %(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    def get_feature_by_topic(self, n_samples=5000) -> pd.DataFrame:
        """
        Gets a Monte-Carlo estimate of the expectation of the feature by topic matrix.

        Parameters
        ----------
        adata
            AnnData to transform. If None, returns the feature by topic matrix for
            the source AnnData.
        n_samples
            Number of samples to take for the Monte-Carlo estimate of the mean.

        Returns
        -------
        A `n_var x n_topics` Pandas DataFrame containing the feature by topic matrix.
        """
        self._check_if_trained(warn=False)

        topic_by_feature = self.module.topic_by_feature(n_samples=n_samples)

        return pd.DataFrame(
            data=topic_by_feature.numpy().T,
            index=self.adata.var_names,
            columns=[f"topic_{i}" for i in range(topic_by_feature.shape[0])],
        )

    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        batch_size: Optional[int] = None,
        n_samples: int = 5000,
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
        n_samples
            Number of samples to take for the Monte-Carlo estimate of the mean.

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
            x = tensors[REGISTRY_KEYS.X_KEY]
            transformed_xs.append(
                self.module.get_topic_distribution(x, n_samples=n_samples)
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
            x = tensors[REGISTRY_KEYS.X_KEY]
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
        total_counts = sum(tensors[REGISTRY_KEYS.X_KEY].sum().item() for tensors in dl)

        return np.exp(
            self.get_elbo(adata=adata, indices=indices, batch_size=batch_size)
            / total_counts
        )
