import logging
from collections.abc import Sequence
from typing import Optional, Union

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from typing_extensions import Literal

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import (
    LayerField,
    ObsmField,
)
from scvi.external.csi.module import Module
from scvi.model.base import BaseModelClass
from scvi.utils import setup_anndata_dsp

from ._training import TrainingCustom
from ._utils import prepare_metadata

logger = logging.getLogger(__name__)


class Model(TrainingCustom, BaseModelClass):
    def __init__(
        self,
        adata: AnnData,
        prior: Literal["standard_normal", "vamp"] = "vamp",
        n_prior_components=5,
        pseudoinputs_data_indices: Optional[np.array] = None,
        **model_kwargs,
    ):
        """CVAE integration model with optional VampPrior and latent cycle-consistency loss

        Parameters
        ----------
        adata
            AnnData object that has been registered via :meth:`~mypackage.MyModel.setup_anndata`.
        prior
            The prior to be used. You can choose between "standard_normal" and "vamp".
        n_prior_components
            Number of prior components in VampPrior.
        pseudoinputs_data_indices
            By default (based on pseudoinputs_data_init),
            VAMP prior pseudoinputs are randomly selected from data.
            Alternatively, one can specify pseudoinput indices using this parameter.
        **model_kwargs
            Keyword args for :class:`~scvi.external.csi.module.Module`
        """
        super().__init__(adata)

        if prior == "vamp":
            if pseudoinputs_data_indices is None:
                pseudoinputs_data_indices = np.random.randint(
                    0, adata.shape[0], n_prior_components
                )
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

        n_cov_const = (
            adata.obsm["covariates"].shape[1] if "covariates" in adata.obsm else 0
        )
        cov_embed_sizes = (
            pd.DataFrame(adata.obsm["covariates_embed"]).nunique(axis=0).to_list()
            if "covariates_embed" in adata.obsm
            else []
        )

        # self.summary_stats provides information about anndata dimensions and other tensor info
        self.module = Module(
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

    @torch.no_grad()
    def embed(
        self,
        adata: AnnData,
        indices: Optional[Sequence[int]] = None,
        cycle: bool = False,
        give_mean: bool = True,
        batch_size: Optional[int] = None,
        as_numpy: bool = True,
    ) -> Union[np.ndarray, torch.Tensor]:
        """
        Return the latent representation for each cell.

        Parameters
        ----------
        adata
            Input adata based on which latent representation is obtained.
        indices
            Data indices to embed. If None embedd all.
        cycle
            Return latent embedding of cycle pass.
        give_mean
            Return posterior mean instead of a sample from posterior.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        as_numpy
            Return iin numpy rather than torch format.

        Returns
        -------
        Latent Embedding
        """
        # Check model and adata
        self._check_if_trained(warn=False)
        adata = self._validate_anndata(adata)
        if indices is None:
            indices = np.arange(adata.n_obs)
        # Prediction
        # Do not shuffle to retain order
        tensors_fwd = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size, shuffle=False
        )
        predicted = []
        for tensors in tensors_fwd:
            # Inference
            inference_inputs = self.module._get_inference_input(tensors)
            inference_outputs = self.module.inference(**inference_inputs)
            if cycle:
                selected_system = self.module.random_select_systems(tensors["system"])
                generative_inputs = self.module._get_generative_input(
                    tensors,
                    inference_outputs,
                    selected_system=selected_system,
                )
                generative_outputs = self.module.generative(
                    **generative_inputs, x_x=False, x_y=True
                )
                inference_cycle_inputs = self.module._get_inference_cycle_input(
                    tensors=tensors,
                    generative_outputs=generative_outputs,
                    selected_system=selected_system,
                )
                inference_outputs = self.module.inference(**inference_cycle_inputs)
            if give_mean:
                predicted += [inference_outputs["z_m"]]
            else:
                predicted += [inference_outputs["z"]]

        predicted = torch.cat(predicted)

        if as_numpy:
            predicted = predicted.cpu().numpy()
        return predicted

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        system_key: str,
        layer: Optional[str] = None,
        categorical_covariate_keys: Optional[list[str]] = None,
        categorical_covariate_embed_keys: Optional[list[str]] = None,
        continuous_covariate_keys: Optional[list[str]] = None,
        covariate_categ_orders: Optional[dict] = None,
        covariate_key_orders: Optional[dict] = None,
        system_order: Optional[list[str]] = None,
        **kwargs,
    ) -> AnnData:
        """
        Prepare adata for input to Model

        Parameters
        ----------
        adata
            Adata object - will be modified in place.
        system_key
            Name of obs column with categorical system information.
        layer
            AnnData layer to use, default X. Should contain normalized and log+1 transformed expression.
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

        # Make sure var names are unique
        if adata.shape[1] != len(set(adata.var_names)):
            raise ValueError("Adata var_names are not unique")

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
            system_order = sorted(adata.obs[system_key].unique())
        # Validate that the provided system_order matches the categories in adata.obs[system_key]
        if set(system_order) != set(adata.obs[system_key].unique()):
            raise ValueError(
                "Provided system_order does not match the categories in adata.obs[system_key]"
            )

        # Make one-hot embedding with specified order
        systems_dict = dict(
            zip(system_order, ([float(i) for i in range(0, len(system_order))]))
        )
        adata.uns["system_order"] = system_order
        system_cat = pd.Series(
            pd.Categorical(
                values=adata.obs[system_key], categories=system_order, ordered=True
            ),
            index=adata.obs.index,
            name="system",
        )
        adata.obsm["system"] = pd.get_dummies(system_cat, dtype=float)

        # Set up covariates
        # TODO this could be handled by specific field type in registry

        # System must not be in cov
        if categorical_covariate_keys is not None:
            if system_key in categorical_covariate_keys:
                raise ValueError("system_key should not be within covariate keys")
        if categorical_covariate_embed_keys is not None:
            if system_key in categorical_covariate_embed_keys:
                raise ValueError("system_key should not be within covariate keys")
        if continuous_covariate_keys is not None:
            if system_key in continuous_covariate_keys:
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
        if (
            continuous_covariate_keys is not None
            or categorical_covariate_keys is not None
        ):
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
        if (
            continuous_covariate_keys is not None
            or categorical_covariate_keys is not None
        ):
            anndata_fields.append(ObsmField("covariates", "covariates"))
        if categorical_covariate_embed_keys is not None:
            anndata_fields.append(ObsmField("covariates_embed", "covariates_embed"))
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
