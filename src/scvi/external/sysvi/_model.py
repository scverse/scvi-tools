from __future__ import annotations

import logging
import warnings
from collections.abc import Sequence
from typing import Literal

import numpy as np
import torch
from anndata import AnnData

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
)
from scvi.model.base import BaseModelClass, UnsupervisedTrainingMixin
from scvi.utils import setup_anndata_dsp

from ._module import SysVAE

logger = logging.getLogger(__name__)


class SysVI(UnsupervisedTrainingMixin, BaseModelClass):
    """Integration model based on cVAE with optional VampPrior and latent cycle-consistency loss.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.external.SysVI.setup_anndata`.
    prior
        The prior distribution to be used. You can choose between "standard_normal" and "vamp".
    n_prior_components
        Number of prior components (i.e. modes) to use in VampPrior.
    pseudoinputs_data_indices
        By default VampPrior pseudoinputs are randomly selected from data.
        Alternatively, one can specify pseudoinput indices using this parameter.
        The number of specified indices in the input 1D array should match n_prior_components
    **model_kwargs
        Keyword args for :class:`~scvi.external.sysvi.SysVAE`
    """

    def __init__(
        self,
        adata: AnnData,
        prior: Literal["standard_normal", "vamp"] = "vamp",
        n_prior_components: int = 5,
        pseudoinputs_data_indices: np.array | None = None,
        **model_kwargs,
    ):
        super().__init__(adata)

        if prior == "vamp":
            if pseudoinputs_data_indices is None:
                pseudoinputs_data_indices = np.random.randint(
                    0, adata.shape[0], n_prior_components
                )
            assert pseudoinputs_data_indices.shape[0] == n_prior_components
            assert pseudoinputs_data_indices.ndim == 1
            pseudoinput_data = next(
                iter(
                    self._make_data_loader(
                        adata=adata,
                        indices=pseudoinputs_data_indices,
                        batch_size=n_prior_components,
                        shuffle=False,
                    )
                )
            )
        else:
            pseudoinput_data = None

        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )

        self.module = SysVAE(
            n_input=self.summary_stats.n_vars,
            n_batch=self.summary_stats.n_batch,
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            n_cats_per_cov=n_cats_per_cov,
            prior=prior,
            n_prior_components=n_prior_components,
            pseudoinput_data=pseudoinput_data,
            **model_kwargs,
        )

        self._model_summary_string = (
            "SysVI - cVAE model with optional VampPrior and latent cycle-consistency loss."
        )
        # necessary line to get params that will be used for saving/loading
        self.init_params_ = self._get_init_params(locals())

        logger.info("The model has been initialized")

    def train(
        self,
        *args,
        plan_kwargs: dict | None = None,
        **kwargs,
    ):
        plan_kwargs = plan_kwargs or {}
        kl_weight_defaults = {"n_epochs_kl_warmup": 0, "n_steps_kl_warmup": 0}
        if any([v != plan_kwargs.get(k, v) for k, v in kl_weight_defaults.items()]):
            warnings.warn(
                "The use of KL weight warmup is not recommended in SysVI. "
                + "The n_epochs_kl_warmup and n_steps_kl_warmup will be reset to 0."
            )
        # Overwrite plan kwargs with kl weight defaults
        plan_kwargs = {**plan_kwargs, **kl_weight_defaults}

        # Pass to parent
        kwargs = kwargs or {}
        kwargs["plan_kwargs"] = plan_kwargs
        super().train(*args, **kwargs)

    @torch.inference_mode()
    def get_latent_representation(
        self,
        adata: AnnData,
        indices: Sequence[int] | None = None,
        give_mean: bool = True,
        batch_size: int | None = None,
        return_dist: bool = False,
    ) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
        """Return the latent representation for each cell.

        Parameters
        ----------
        adata
            Input adata for which latent representation should be obtained.
        indices
            Data indices to embed. If None embedd all cells.
        give_mean
            Return the posterior mean instead of a sample from the posterior.
            Ignored if `return_dist` is ``True``.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_dist
            If ``True``, returns the mean and variance of the latent distribution. Otherwise,
            returns the mean of the latent distribution.

        Returns
        -------
        Latent Embedding
        """
        # Check model and adata
        self._check_if_trained(warn=False)
        adata = self._validate_anndata(adata)
        if indices is None:
            indices = np.arange(adata.n_obs)
        # Do not shuffle to retain order
        tensors_fwd = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size, shuffle=False
        )
        predicted_m = []
        predicted_v = []
        for tensors in tensors_fwd:
            # Inference
            inference_inputs = self.module._get_inference_input(tensors)
            inference_outputs = self.module.inference(**inference_inputs)
            if give_mean or return_dist:
                predicted_m += [inference_outputs["z_m"]]
            else:
                predicted_m += [inference_outputs["z"]]
            if return_dist:
                predicted_v += [inference_outputs["z_v"]]

        predicted_m = torch.cat(predicted_m).cpu().numpy()
        if return_dist:
            predicted_v = torch.cat(predicted_v).cpu().numpy()

        if return_dist:
            return predicted_m, predicted_v
        else:
            return predicted_m

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        batch_key: str,
        layer: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        weight_batches: bool = False,
        **kwargs,
    ):
        """Prepare adata for input to Model

        Setup distinguishes between two main types of covariates that can be corrected for:
        - batch (referred to as "batch" in the original publication): Single categorical covariate that
        will be corrected via cycle consistency loss. It will be also used as a condition in cVAE.
        This covariate is expected to correspond to stronger batch effects, such as between datasets from the
        different sequencing technology or model systems (animal species, in-vitro models, etc.).
        - covariate (includes both continous and categorical covariates): Additional covariates to be used only
        as a condition in cVAE, but not corrected via cycle loss.
        These covariates are expected to correspond to weaker batch effects, such as between datasets from the
        same sequencing technology and model batch (animal, in-vitro, etc.) or between samples within a dataset.

        Parameters
        ----------
        adata
            Adata object - will be modified in place.
        batch_key
            Name of the obs column with the substantial batch effect covariate,
            referred to as batch in the original publication (Hrovatin, et al., 2023).
            Should be categorical.
        layer
            AnnData layer to use, default is X.
            Should contain normalized and log+1 transformed expression.
        categorical_covariate_keys
            Name of obs column with additional categorical covariate information.
            Will be one hot encoded or embedded, as later defined in the model.
        continuous_covariate_keys
            Name of obs column with additional continuous covariate information.
        """
        setup_method_args = cls._get_setup_method_args(**locals())

        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=False),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]
        if weight_batches:
            warnings.warn("The use of inverse batch proportion weights is experimental.")
            batch_weights_key = "batch_weights"
            adata.obs[batch_weights_key] = adata.obs[batch_key].map(
                {cat: 1 / n for cat, n in adata.obs[batch_key].value_counts().items()}
            )
            anndata_fields.append(NumericalObsField(batch_weights_key, batch_weights_key))
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
