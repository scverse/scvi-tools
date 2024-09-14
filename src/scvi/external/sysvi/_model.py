from __future__ import annotations

import logging
from collections.abc import Sequence
from typing import Literal, Tuple

import numpy as np
import pandas as pd
import torch
from anndata import AnnData

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data._constants import _SCVI_UUID_KEY
from scvi.data._utils import _check_if_view
from scvi.data.fields import (
    LayerField,
    ObsmField,
)
from scvi.model.base import BaseModelClass, UnsupervisedTrainingMixin
from scvi.utils import setup_anndata_dsp

from ._module import SysVAE
from ._trainingplans import TrainingPlanCustom
from ._utils import prepare_metadata

logger = logging.getLogger(__name__)


class TrainingCustom(UnsupervisedTrainingMixin):
    """Train method with custom TrainingPlan."""

    # TODO could make custom Trainer (in a custom TrainRunner) to have in init params for early stopping
    #  "loss" rather than "elbo" components in available param specifications - for now just use
    #  a loss that is against the param specification

    _training_plan_cls = TrainingPlanCustom


class SysVI(TrainingCustom, BaseModelClass):
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

        n_cov_const = adata.obsm["covariates"].shape[1] if "covariates" in adata.obsm else 0
        cov_embed_sizes = (
            pd.DataFrame(adata.obsm["covariates_embed"]).nunique(axis=0).to_list()
            if "covariates_embed" in adata.obsm
            else []
        )

        self.module = SysVAE(
            n_input=adata.shape[1],
            n_cov_const=n_cov_const,
            cov_embed_sizes=cov_embed_sizes,
            n_system=adata.obsm["system"].shape[1],
            prior=prior,
            n_prior_components=n_prior_components,
            pseudoinput_data=pseudoinput_data,
            **model_kwargs,
        )

        self._model_summary_string = (
            "cVAE model with optional VampPrior and latent cycle-consistency loss"
        )
        # necessary line to get params that will be used for saving/loading
        self.init_params_ = self._get_init_params(locals())

        logger.info("The model has been initialized")

    @torch.inference_mode()
    def get_latent_representation(
        self,
        adata: AnnData,
        indices: Sequence[int] | None = None,
        give_mean: bool = True,
        batch_size: int | None = None,
        return_dist: bool = False,
    ) -> np.ndarray | Tuple[np.ndarray, np.ndarray]:
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

    def _validate_anndata(
        self, adata: AnnData | None = None, copy_if_view: bool = True
    ) -> AnnData:
        """Validate anndata has been properly registered

        Parameters
        ----------
        adata
            Adata to validate. If None use SysVI's adata.
        copy_if_view
            Whether to copy adata before

        Returns
        -------

        """
        if adata is None:
            adata = self.adata

        _check_if_view(adata, copy_if_view=copy_if_view)

        if _SCVI_UUID_KEY not in adata.uns:
            raise ValueError("Adata is not set up. Use SysVI.setup_anndata first.")
        else:
            # Check that all required fields are present and match the Model's adata
            assert (
                self.adata.uns["layer_information"]["layer"]
                == adata.uns["layer_information"]["layer"]
            )
            assert (
                self.adata.uns["layer_information"]["var_names"]
                == adata.uns["layer_information"]["var_names"]
            )
            assert self.adata.uns["system_order"] == adata.uns["system_order"]
            for covariate_type, covariate_keys in self.adata.uns["covariate_key_orders"].items():
                assert covariate_keys == adata.uns["covariate_key_orders"][covariate_type]
                if "categorical" in covariate_type:
                    for covariate_key in covariate_keys:
                        assert (
                            self.adata.uns["covariate_categ_orders"][covariate_key]
                            == adata.uns["covariate_categ_orders"][covariate_key]
                        )
        # Ensures that manager is set up
        super()._validate_anndata(adata)

        return adata

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        batch_key: str,
        layer: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        categorical_covariate_embed_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        covariate_categ_orders: dict | None = None,
        covariate_key_orders: dict | None = None,
        system_order: list[str] | None = None,
        **kwargs,
    ):
        """Prepare adata for input to Model

        Parameters
        ----------
        adata
            Adata object - will be modified in place.
        batch_key
            Name of the obs column with the substantial batch effect covariate,
            referred to as system in the original publication (Hrovatin, et al., 2023).
            Should be categorical.
        layer
            AnnData layer to use, default is X.
            Should contain normalized and log+1 transformed expression.
        categorical_covariate_keys
            Name of obs column with additional categorical covariate information. Will be one hot encoded.
        categorical_covariate_embed_keys
            Name of obs column with additional categorical covariate information. Embedding will be learned.
            This can be useful if the number of categories is very large, which would increase memory usage.
            If using this type of covariate representation please also cite
            `scPoli <[https://doi.org/10.1038/s41592-023-02035-2]>`_ .
        continuous_covariate_keys
            Name of obs column with additional continuous covariate information.
        covariate_categ_orders
            Covariate encoding information. Should be used if a new adata is to be set up according
            to setup of an existing adata. Access via adata.uns['covariate_categ_orders'] of already setup adata.
        covariate_key_orders
            Covariate encoding information. Should be used if a new adata is to be set up according
            to setup of an existing adata. Access via adata.uns['covariate_key_orders'] of already setup adata.
        system_order
            Same as covariate_orders, but for system. Access via adata.uns['system_order']
        """
        setup_method_args = cls._get_setup_method_args(**locals())

        if adata.shape[1] != len(set(adata.var_names)):
            raise ValueError("Adata var_names are not unique")

        # The used layer argument
        # This could be also done via registry, but that is too cumbersome
        adata.uns["layer_information"] = {
            "layer": layer,
            "var_names": list(adata.var_names),
        }

        # If setup is to be prepared wtr another adata specs make sure all relevant info is present
        if covariate_categ_orders or covariate_key_orders or system_order:
            assert system_order is not None
            if (
                categorical_covariate_keys is not None
                or categorical_covariate_embed_keys is not None
                or continuous_covariate_keys is not None
            ):
                assert covariate_categ_orders is not None
                assert covariate_key_orders is not None

        # Make system embedding with specific category order

        # Define order of system categories
        if system_order is None:
            system_order = sorted(adata.obs[batch_key].unique())
        # Validate that the provided system_order matches the categories in adata.obs[system_key]
        if set(system_order) != set(adata.obs[batch_key].unique()):
            raise ValueError(
                "Provided system_order does not match the categories in adata.obs[system_key]"
            )

        # Make one-hot embedding with specified order
        systems_dict = dict(zip(system_order, ([float(i) for i in range(0, len(system_order))])))
        adata.uns["system_order"] = system_order
        system_cat = pd.Series(
            pd.Categorical(values=adata.obs[batch_key], categories=system_order, ordered=True),
            index=adata.obs.index,
            name="system",
        )
        adata.obsm["system"] = pd.get_dummies(system_cat, dtype=float)

        # Set up covariates
        # TODO this could be handled by specific field type in registry

        # System must not be in cov
        if categorical_covariate_keys is not None:
            if batch_key in categorical_covariate_keys:
                raise ValueError("system_key should not be within covariate keys")
        if categorical_covariate_embed_keys is not None:
            if batch_key in categorical_covariate_embed_keys:
                raise ValueError("system_key should not be within covariate keys")
        if continuous_covariate_keys is not None:
            if batch_key in continuous_covariate_keys:
                raise ValueError("system_key should not be within covariate keys")

        # Prepare covariate training representations/embedding
        covariates, covariates_embed, orders_dict, cov_dict = prepare_metadata(
            meta_data=adata.obs,
            cov_cat_keys=categorical_covariate_keys,
            cov_cat_embed_keys=categorical_covariate_embed_keys,
            cov_cont_keys=continuous_covariate_keys,
            categ_orders=covariate_categ_orders,
            key_orders=covariate_key_orders,
        )

        # Save covariate representation and order information
        adata.uns["covariate_categ_orders"] = orders_dict
        adata.uns["covariate_key_orders"] = cov_dict
        if continuous_covariate_keys is not None or categorical_covariate_keys is not None:
            adata.obsm["covariates"] = covariates
        else:
            # Remove if present since the presence of this key
            # is in model used to determine if cov should be used or not
            if "covariates" in adata.obsm:
                del adata.obsm["covariates"]
        if categorical_covariate_embed_keys is not None:
            adata.obsm["covariates_embed"] = covariates_embed
        else:
            if "covariates_embed" in adata.obsm:
                del adata.obsm["covariates_embed"]

        # Anndata setup

        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=False),
            ObsmField("system", "system"),
        ]
        # Covariate fields are optional
        if continuous_covariate_keys is not None or categorical_covariate_keys is not None:
            anndata_fields.append(ObsmField("covariates", "covariates"))
        if categorical_covariate_embed_keys is not None:
            anndata_fields.append(ObsmField("covariates_embed", "covariates_embed"))
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
