from __future__ import annotations

import logging
import warnings
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData

import numpy as np

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
)
from scvi.model.base import BaseModelClass, RNASeqMixin, UnsupervisedTrainingMixin, VAEMixin
from scvi.utils import setup_anndata_dsp

from ._module import SysVAE

logger = logging.getLogger(__name__)


class SysVI(UnsupervisedTrainingMixin, RNASeqMixin, VAEMixin, BaseModelClass):
    """Integration with cVAE & optional VampPrior and latent cycle-consistency.

     Described in
     `Hrovatin et al. (2023) <https://doi.org/10.1101/2023.11.03.565463>`_.

    Parameters
    ----------
    adata
        AnnData object that has been registered via
        :meth:`~scvi.external.SysVI.setup_anndata`.
    prior
        The prior distribution to be used.
        You can choose between ``"standard_normal"`` and ``"vamp"``.
    n_prior_components
        Number of prior components (i.e. modes) to use in VampPrior.
    pseudoinputs_data_indices
        By default, VampPrior pseudoinputs are randomly selected from data.
        Alternatively, one can specify pseudoinput indices using this parameter.
        The number of specified indices in the input 1D array should match
        ``n_prior_components``.
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
                    0, self.summary_stats.n_vars, n_prior_components
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
            "SysVI - cVAE model with optional VampPrior " + "and latent cycle-consistency loss."
        )
        self.init_params_ = self._get_init_params(locals())

        logger.info("The model has been initialized")

    def train(
        self,
        plan_kwargs: dict | None = None,
        **train_kwargs,
    ):
        """Train the models.

        Overwrites the ``train`` method of
        class:`~scvi.model.base.UnsupervisedTrainingMixin`
        to prevent the use of KL loss warmup (specified in ``plan_kwargs``).
        This is disabled as our experiments showed poor integration in the
        cycle model when using KL loss warmup.

        Parameters
        ----------
        plan_kwargs
            Training plan kwargs in `meth`:`scvi.train.TrainingPlan`.
        train_kwargs
            Training kwargs. Passed to `meth`:`scvi.model.base.BaseModelClass.train`.
        """
        plan_kwargs = plan_kwargs or {}
        kl_weight_defaults = {"n_epochs_kl_warmup": 0, "n_steps_kl_warmup": 0}
        if any(v != plan_kwargs.get(k, v) for k, v in kl_weight_defaults.items()):
            warnings.warn(
                "The use of KL weight warmup is not recommended in SysVI. "
                + "The n_epochs_kl_warmup and n_steps_kl_warmup "
                + "will be reset to 0.",
                stacklevel=settings.warnings_stacklevel,
            )
        plan_kwargs["n_epochs_kl_warmup"] = 0
        plan_kwargs["n_steps_kl_warmup"] = 0

        # Pass to parent
        train_kwargs = train_kwargs or {}
        train_kwargs["plan_kwargs"] = plan_kwargs
        super().train(**train_kwargs)

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
        """Prepare adata for input to SysVI model.

        Setup distinguishes between two main types of covariates that can be
        corrected for:

        - batch - referred to as "system" in the original publication
          Hrovatin, et al., 2023):
          Single categorical covariate that will be corrected via cycle
          consistency loss.
          It will be also used as a condition in cVAE.
          This covariate is expected to correspond to stronger batch effects,
          such as between datasets from different sequencing technology or
          model systems (animal species, in-vitro models and tissue, etc.).
        - covariate (includes both continous and categorical covariates):
          Additional covariates to be used only
          as a condition in cVAE, but not corrected via cycle loss.
          These covariates are expected to correspond to weaker batch effects,
          such as between datasets from the same sequencing technology and
          system (animal, in-vitro, etc.) or between samples within a dataset.

        Parameters
        ----------
        %(param_adata)s
        %(param_batch_key)s
        layer
            AnnData layer to use, default is X.
            Should contain normalized and log1p transformed expression.
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())

        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=False),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]
        if weight_batches:
            warnings.warn(
                "The use of inverse batch proportion weights " + "is experimental.",
                stacklevel=settings.warnings_stacklevel,
            )
            batch_weights_key = "batch_weights"
            adata.obs[batch_weights_key] = adata.obs[batch_key].map(
                {cat: 1 / n for cat, n in adata.obs[batch_key].value_counts().items()}
            )
            anndata_fields.append(NumericalObsField(batch_weights_key, batch_weights_key))
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
