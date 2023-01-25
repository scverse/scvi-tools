import logging
from typing import Literal, Optional, Union

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from torch.distributions.multinomial import Multinomial

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalObsField, LayerField, NumericalObsField
from scvi.model._utils import _init_library_size
from scvi.model.base import (
    BaseModelClass,
    RNASeqMixin,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.utils import setup_anndata_dsp, track

from ._module import SCAR_VAE

logger = logging.getLogger(__name__)


class SCAR(RNASeqMixin, VAEMixin, UnsupervisedTrainingMixin, BaseModelClass):
    """
    Ambient RNA removal in scRNA-seq data :cite:p:`Sheng22`.

    Original Github: https://github.com/Novartis/scar.
    The models are parameter matched in architecture, activations, dropout, sparsity, and batch normalization.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi_external.SCAR.setup_anndata`.
    ambient_profile
        The probability of occurrence of each ambient transcript.\
            If None, averaging cells to estimate the ambient profile, by default None.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    gene_likelihood
        One of:
        * ``'b'`` - Binomial distribution
        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution
    latent_distribution
        One of:
        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
    scale_activation
        Activation layer to use for px_scale_decoder
    sparsity
        The sparsity of expected native signals. It varies between datasets,
        e.g. if one prefilters genes -- use only highly variable genes --
        the sparsity should be low; on the other hand, it should be set high
        in the case of unflitered genes.
    **model_kwargs
        Keyword args for :class:`~scvi_external.SCAR`
    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> raw_adata = anndata.read_h5ad(path_to_raw_anndata)
    >>> scvi_external.SCAR.setup_anndata(adata, batch_key="batch")
    >>> scvi_external.SCAR.get_ambient_profile(adata=adata, raw_adata=raw_adata, prob=0.995)
    >>> vae = scvi_external.SCAR(adata)
    >>> vae.train()
    >>> adata.obsm["X_scAR"] = vae.get_latent_representation()
    >>> adata.layers['denoised'] = vae.get_denoised_counts()
    """

    def __init__(
        self,
        adata: AnnData,
        ambient_profile: Union[str, np.ndarray, pd.DataFrame, torch.tensor] = None,
        n_hidden: int = 150,
        n_latent: int = 15,
        n_layers: int = 2,
        dropout_rate: float = 0.0,
        gene_likelihood: Literal["zinb", "nb", "b", "poisson"] = "b",
        latent_distribution: Literal["normal", "ln"] = "normal",
        scale_activation: Literal["softmax", "softplus", "softplus_sp"] = "softplus_sp",
        sparsity: float = 0.9,
        **model_kwargs,
    ):
        super().__init__(adata)

        n_batch = self.summary_stats.n_batch
        use_size_factor_key = (
            REGISTRY_KEYS.SIZE_FACTOR_KEY in self.adata_manager.data_registry
        )
        library_log_means, library_log_vars = None, None
        if not use_size_factor_key:
            library_log_means, library_log_vars = _init_library_size(
                self.adata_manager, n_batch
            )

        # self.summary_stats provides information about anndata dimensions and other tensor info
        if not torch.is_tensor(ambient_profile):
            if isinstance(ambient_profile, str):
                ambient_profile = np.nan_to_num(adata.varm[ambient_profile])
            elif isinstance(ambient_profile, pd.DataFrame):
                ambient_profile = ambient_profile.fillna(0).values
            elif isinstance(ambient_profile, np.ndarray):
                ambient_profile = np.nan_to_num(ambient_profile)
            elif not ambient_profile:
                ambient_profile = adata.X.sum(axis=0) / adata.X.sum(axis=0).sum()
                ambient_profile = np.nan_to_num(ambient_profile)
            else:
                raise TypeError(
                    f"Expecting str / np.array / None / pd.DataFrame, but get a {type(ambient_profile)}"
                )
            ambient_profile = (
                torch.from_numpy(np.asarray(ambient_profile)).float().reshape(1, -1)
            )

        self.module = SCAR_VAE(
            ambient_profile=ambient_profile,
            n_input=self.summary_stats.n_vars,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
            use_size_factor_key=use_size_factor_key,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            scale_activation=scale_activation,
            sparsity=sparsity,
            **model_kwargs,
        )
        self._model_summary_string = (
            "SCVI-AR Model with the following params: \nn_hidden: {}, n_latent: {}, n_layers: {}, dropout_rate: "
            "{}, gene_likelihood: {}, latent_distribution: {}"
        ).format(
            n_hidden,
            n_latent,
            n_layers,
            dropout_rate,
            gene_likelihood,
            latent_distribution,
        )
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: Optional[str] = None,
        size_factor_key: Optional[str] = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
        %(param_layer)s
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_size_factor_key)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, None),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, None),
            NumericalObsField(
                REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False
            ),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    @staticmethod
    def get_ambient_profile(
        adata: AnnData,
        raw_adata: AnnData,
        prob: float = 0.995,
        min_raw_counts: int = 2,
        iterations: int = 3,
        n_batch: int = 1,
        sample: int = 50000,
    ):
        """
        Calculate ambient profile for relevant features.

        Identify the cell-free droplets through a multinomial distribution. See EmptyDrops :cite:p:`Lun2019` for details.

        Parameters
        ----------
        adata
            A filtered adata object, loaded from filtered_feature_bc_matrix using `scanpy.read`, gene filtering is
            recommended to save memory.
        raw_adata
            A raw adata object, loaded from raw_feature_bc_matrix using :meth:`~scanpy.read`.
        prob
            The probability of each gene, considered as containing ambient RNA if greater than prob
            (joint prob euqals to the product of all genes for a droplet), by default 0.995.
        min_raw_counts
            Total counts filter for raw_adata, filtering out low counts to save memory, by default 2.
        iterations
            Total iterations, by default 3.
        n_batch
            Total number of batches, set it to a bigger number when out of memory issue occurs, by default 1.
        sample
            Randomly sample droplets to test, if greater than total droplets, use all droplets, by default 50000.

        Returns
        -------
        The relevant ambient profile is added in `adata.varm`
        """
        # take subset genes to save memory
        raw_adata = raw_adata[:, raw_adata.var_names.isin(adata.var_names)]
        raw_adata = raw_adata[raw_adata.X.sum(axis=1) >= min_raw_counts]

        sample = int(sample)
        idx = np.random.choice(
            raw_adata.shape[0], size=min(raw_adata.shape[0], sample), replace=False
        )
        raw_adata = raw_adata[idx]
        print(f"Randomly sampling {sample} droplets to calculate the ambient profile.")
        # initial estimation of ambient profile, will be updated
        ambient_prof = raw_adata.X.sum(axis=0) / raw_adata.X.sum()

        for _ in track(range(iterations)):
            # calculate joint probability (log) of being cell-free droplets for each droplet
            log_prob = []
            batch_idx = np.floor(
                np.array(range(raw_adata.shape[0])) / raw_adata.shape[0] * n_batch
            )
            for b in range(n_batch):
                try:
                    count_batch = raw_adata[batch_idx == b].X.astype(int).A
                except MemoryError:
                    raise MemoryError("use more batches by setting a higher n_batch")
                log_prob_batch = Multinomial(
                    probs=torch.tensor(ambient_prof), validate_args=False
                ).log_prob(torch.Tensor(count_batch))
                log_prob.append(log_prob_batch)
            log_prob = np.concatenate(log_prob, axis=0)
            raw_adata.obs["log_prob"] = log_prob
            raw_adata.obs["droplets"] = "other droplets"
            # cell-containing droplets
            raw_adata.obs.loc[
                raw_adata.obs_names.isin(adata.obs_names), "droplets"
            ] = "cells"
            # identify cell-free droplets
            raw_adata.obs["droplets"] = raw_adata.obs["droplets"].mask(
                raw_adata.obs["log_prob"] >= np.log(prob) * raw_adata.shape[1],
                "cell-free droplets",
            )
            emptydrops = raw_adata[raw_adata.obs["droplets"] == "cell-free droplets"]
            if emptydrops.shape[0] < 50:
                raise Exception("Too few emptydroplets. Lower the prob parameter")
            ambient_prof = emptydrops.X.sum(axis=0) / emptydrops.X.sum()

        # update ambient profile
        adata.varm["ambient_profile"] = np.asarray(
            emptydrops.X.sum(axis=0).reshape(-1, 1) / emptydrops.X.sum()
        )

    @torch.no_grad()
    def get_denoised_counts(
        self,
        adata: Optional[AnnData] = None,
        n_samples: int = 1,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        r"""
        Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        n_samples
            Number of samples for each cell.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        x_denoised : :py:class:`torch.Tensor`
            tensor with shape (n_cells, n_genes)
        """
        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(adata=adata, batch_size=batch_size)

        data_loader_list = []
        for tensors in scdl:
            x = tensors[REGISTRY_KEYS.X_KEY]
            inference_kwargs = dict(n_samples=n_samples)
            _, generative_outputs = self.module.forward(
                tensors=tensors,
                inference_kwargs=inference_kwargs,
                compute_loss=False,
            )

            total_count_per_cell = x.sum(dim=1).reshape(-1, 1)

            if self.module.gene_likelihood == "b":
                px_scale = generative_outputs["px"].probs
            else:
                px_scale = generative_outputs["px"].scale
            expected_counts = total_count_per_cell * px_scale.cpu()

            b = torch.distributions.Binomial(
                probs=expected_counts - expected_counts.floor()
            )
            expected_counts = expected_counts.floor() + b.sample()

            if n_samples > 1:
                expected_counts = torch.median(expected_counts, dim=0)[0]

            data_loader_list.append(expected_counts)

        x_denoised = torch.cat(data_loader_list, dim=0).numpy()

        return x_denoised
