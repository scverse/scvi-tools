from __future__ import annotations

import logging
from collections import defaultdict
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable

    from anndata import AnnData
    from mudata import MuData

import numpy as np
import sparse
import torch

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager, fields
from scvi.data._constants import _SETUP_ARGS_KEY
from scvi.external.methylvi._base_components import BSSeqMixin
from scvi.external.methylvi._utils import _context_cov_key, _context_mc_key
from scvi.model.base import (
    ArchesMixin,
    BaseModelClass,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.utils import setup_anndata_dsp

from ._methylvi_module import METHYLVAE

logger = logging.getLogger(__name__)


class METHYLVI(VAEMixin, BSSeqMixin, UnsupervisedTrainingMixin, ArchesMixin, BaseModelClass):
    """
    Model class for methylVI :cite:p:`Weinberger2023a`

    Parameters
    ----------
    mdata
        MuData object that has been registered via :meth:`~scvi.external.METHYLVI.setup_mudata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    **model_kwargs
        Keyword args for :class:`~scvi.external.methylvi.METHYLVAE`

    Examples
    --------
    >>> mdata = mudata.read_h5mu(path_to_mudata)
    >>> MethylVI.setup_mudata(mdata, batch_key="batch")
    >>> vae = MethylVI(mdata)
    >>> vae.train()
    >>> mdata.obsm["X_methylVI"] = vae.get_latent_representation()
    """

    def __init__(
        self,
        mdata: MuData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        **model_kwargs,
    ):
        super().__init__(mdata)

        n_batch = self.summary_stats.n_batch
        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY)[
                fields.CategoricalJointObsField.N_CATS_PER_KEY
            ]
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )

        self.contexts = self.get_anndata_manager(mdata, required=True).registry[_SETUP_ARGS_KEY][
            "methylation_contexts"
        ]
        self.num_features_per_context = [mdata[context].shape[1] for context in self.contexts]
        self.get_normalized_function_name = "get_normalized_methylation"

        n_input = np.sum(self.num_features_per_context)

        self.module = METHYLVAE(
            n_input=n_input,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            n_batch=n_batch,
            n_cats_per_cov=n_cats_per_cov,
            contexts=self.contexts,
            num_features_per_context=self.num_features_per_context,
            **model_kwargs,
        )
        self._model_summary_string = (
            "Overwrite this attribute to get an informative representation for your model"
        )
        # necessary line to get params that will be used for saving/loading
        self.init_params_ = self._get_init_params(locals())

        logger.info("The model has been initialized")

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        **kwargs,
    ) -> AnnData | None:
        """
        %(summary)s.

        Parameters
        ----------
        %(param_adata)s

        Returns
        -------
        %(returns)s
        """
        raise NotImplementedError("METHYLVI must be used with a MuData object.")

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_mudata(
        cls,
        mdata: MuData,
        mc_layer: str,
        cov_layer: str,
        methylation_contexts: Iterable[str],
        batch_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        modalities=None,
        **kwargs,
    ):
        """%(summary_mdata)s.

        Parameters
        ----------
        %(param_mdata)s
        mc_layer
            Layer containing methylated cytosine counts for each set of methylation features.
        cov_layer
            Layer containing total coverage counts for each set of methylation features.
        methylation_contexts
            List of modality fields in `mdata` object representing different methylation contexts.
            Each context must be equipped with a layer containing the number of methylated counts
            (specified by `mc_layer`) and total number of counts (specified by `cov_layer`) for
            each genomic region feature.
        %(param_batch_key)s
        %(param_cat_cov_keys)s
        %(param_modalities)s

        Examples
        --------
        MethylVI.setup_mudata(
            mdata,
            mc_layer="mc",
            cov_layer="cov",
            batch_key="Platform",
            methylation_modalities=['mCG', 'mCH'],
            modalities={
                "batch_key": "mCG"
            },
        )

        """
        if modalities is None:
            modalities = {}
        setup_method_args = METHYLVI._get_setup_method_args(**locals())

        if methylation_contexts is None:
            raise ValueError("Methylation contexts cannot be None.")

        modalities_ = cls._create_modalities_attr_dict(modalities, setup_method_args)

        batch_field = fields.MuDataCategoricalObsField(
            REGISTRY_KEYS.BATCH_KEY,
            batch_key,
            mod_key=modalities_.batch_key,
        )

        cat_cov_field = fields.MuDataCategoricalJointObsField(
            REGISTRY_KEYS.CAT_COVS_KEY,
            categorical_covariate_keys,
            mod_key=modalities_.categorical_covariate_keys,
        )

        mc_fields = []
        cov_fields = []

        for context in methylation_contexts:
            mc_fields.append(
                fields.MuDataLayerField(
                    _context_mc_key(context),
                    mc_layer,
                    mod_key=context,
                    is_count_data=True,
                    mod_required=True,
                )
            )

            cov_fields.append(
                fields.MuDataLayerField(
                    _context_cov_key(context),
                    cov_layer,
                    mod_key=context,
                    is_count_data=True,
                    mod_required=True,
                )
            )

        mudata_fields = mc_fields + cov_fields + [batch_field] + [cat_cov_field]
        adata_manager = AnnDataManager(fields=mudata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(mdata, **kwargs)

        cls.register_manager(adata_manager)

    @torch.inference_mode()
    def posterior_predictive_sample(
        self,
        mdata: MuData | None = None,
        n_samples: int = 1,
        batch_size: int | None = None,
    ) -> dict[str, sparse.GCXS] | sparse.GCXS:
        r"""
        Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        mdata
            MuData object with equivalent structure to initial MuData. If `None`, defaults to the
            MuData object used to initialize the model.
        n_samples
            Number of samples for each cell.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        x_new : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_regions, n_samples)
        """
        mdata = self._validate_anndata(mdata)

        scdl = self._make_data_loader(adata=mdata, batch_size=batch_size)

        x_new = defaultdict(list)
        for tensors in scdl:
            samples = self.module.sample(
                tensors,
                n_samples=n_samples,
            )

            for context in self.contexts:
                x_new[context].append(sparse.GCXS.from_numpy(samples[context].numpy()))

        for context in self.contexts:
            x_new[context] = sparse.concatenate(
                x_new[context]
            )  # Shape (n_cells, n_regions, n_samples)

        return x_new
