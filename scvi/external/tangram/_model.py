import logging
from typing import Dict, Literal, Optional, Union

import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd
import scipy
from anndata import AnnData
from jaxlib.xla_extension import Device
from mudata import MuData

from scvi.data import AnnDataManager, AnnDataManagerValidationCheck, fields
from scvi.external.tangram._module import TANGRAM_REGISTRY_KEYS, TangramMapper
from scvi.model.base import BaseModelClass
from scvi.train import JaxTrainingPlan
from scvi.utils import setup_anndata_dsp, track

logger = logging.getLogger(__name__)


def _asarray(x: np.ndarray, device: Device) -> jnp.ndarray:
    return jax.device_put(x, device=device)


class Tangram(BaseModelClass):
    """
    Reimplementation of Tangram :cite:p:`Biancalani21` for mapping single-cell RNA-seq data to spatial data.

    Currently the "cells" and "constrained" modes are implemented.

    Original code:
    https://github.com/broadinstitute/Tangram.

    Parameters
    ----------
    mdata
        MuData object that has been registered via :meth:`~scvi.external.Tangram.setup_mudata`.
    constrained
        Whether to use the constrained version of Tangram instead of cells mode.
    target_count
        The number of cells to be filtered. Necessary when `constrained` is True.
    **model_kwargs
        Keyword args for :class:`~scvi.external.tangram.TangramMapper`

    Examples
    --------
    >>> from scvi.external import Tangram
    >>> ad_sc = anndata.read_h5ad(path_to_sc_anndata)
    >>> ad_sp = anndata.read_h5ad(path_to_sp_anndata)
    >>> markers = pd.read_csv(path_to_markers, index_col=0) # genes to use for mapping
    >>> mdata = mudata.MuData({"sp_full": ad_sp, "sc_full": ad_sc, "sp": ad_sp[:, markers].copy(), "sc": ad_sc[:, markers].copy()})
    >>> modalities = {"density_prior_key": "sp", "sc_layer": "sc", "sp_layer": "sp"}
    >>> Tangram.setup_mudata(mdata, density_prior_key="rna_count_based_density", modalities=modalities)
    >>> tangram = Tangram(sc_adata)
    >>> tangram.train()
    >>> ad_sc.obsm["tangram_mapper"] = tangram.get_mapper_matrix()
    >>> ad_sp.obsm["tangram_cts"] = tangram.project_cell_annotations(ad_sc, ad_sp, ad_sc.obsm["tangram_mapper"], ad_sc.obs["labels"])
    >>> projected_ad_sp = tangram.project_genes(ad_sc, ad_sp, ad_sc.obsm["tangram_mapper"])

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/tangram_scvi_tools`
    """

    def __init__(
        self,
        sc_adata: AnnData,
        constrained: bool = False,
        target_count: Optional[int] = None,
        **model_kwargs,
    ):
        super().__init__(sc_adata)
        self.n_obs_sc = self.adata_manager.get_from_registry(
            TANGRAM_REGISTRY_KEYS.SC_KEY
        ).shape[0]
        self.n_obs_sp = self.adata_manager.get_from_registry(
            TANGRAM_REGISTRY_KEYS.SP_KEY
        ).shape[0]

        if constrained and target_count is None:
            raise ValueError(
                "Please specify `target_count` when using constrained Tangram."
            )
        has_density_prior = not self.adata_manager.fields[-1].is_empty
        if has_density_prior:
            prior = self.adata_manager.get_from_registry(
                TANGRAM_REGISTRY_KEYS.DENSITY_KEY
            )
            if np.abs(prior.ravel().sum() - 1) > 1e-3:
                raise ValueError(
                    "Density prior must sum to 1. Please normalize the prior."
                )

        self.module = TangramMapper(
            n_obs_sc=self.n_obs_sc,
            n_obs_sp=self.n_obs_sp,
            lambda_d=1.0 if has_density_prior else 0.0,
            constrained=constrained,
            target_count=target_count,
            **model_kwargs,
        )
        self._model_summary_string = (
            "TangramMapper Model with params: \nn_obs_sc: {}, n_obs_sp: {}"
        ).format(
            self.n_obs_sc,
            self.n_obs_sp,
        )
        self.init_params_ = self._get_init_params(locals())

    def get_mapper_matrix(self) -> np.ndarray:
        """
        Return the mapping matrix.

        Returns
        -------
        Mapping matrix of shape (n_obs_sp, n_obs_sc)
        """
        return jax.device_get(
            jax.nn.softmax(self.module.params["mapper_unconstrained"], axis=1)
        )

    def train(
        self,
        max_epochs: int = 1000,
        use_gpu: Optional[Union[str, int, bool]] = None,
        lr: float = 0.1,
        plan_kwargs: Optional[dict] = None,
    ):
        """
        Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset.
        use_gpu
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
        lr
            Optimiser learning rate (default optimiser is :class:`~pyro.optim.ClippedAdam`).
            Specifying optimiser via plan_kwargs overrides this choice of lr.
        plan_kwargs
            Keyword args for :class:`~scvi.train.JaxTrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        """
        update_dict = {
            "optim_kwargs": {
                "learning_rate": lr,
                "eps": 1e-8,
                "weight_decay": 0,
            }
        }
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict
        device = jax.devices("cpu")[0]
        if use_gpu is None or use_gpu is True:
            try:
                device = jax.devices("gpu")[0]
                self.module.to(device)
                logger.info(
                    "Jax module moved to GPU. "
                    "Note: Pytorch lightning will show GPU is not being used for the Trainer."
                )
            except RuntimeError:
                logger.debug("No GPU available to Jax.")
        else:
            self.module.to(device)
            logger.info("Jax module moved to CPU.")
        tensor_dict = self._get_tensor_dict(device=device)
        training_plan = JaxTrainingPlan(self.module, **plan_kwargs)
        module_init = self.module.init(self.module.rngs, tensor_dict)
        state, params = module_init.pop("params")
        training_plan.set_train_state(params, state)
        train_step_fn = JaxTrainingPlan.jit_training_step
        pbar = track(range(max_epochs), style="tqdm", description="Training")
        history = pd.DataFrame(index=np.arange(max_epochs), columns=["loss"])
        for i in pbar:
            self.module.train_state, loss, _ = train_step_fn(
                self.module.train_state, tensor_dict, self.module.rngs
            )
            loss = jax.device_get(loss)
            history.iloc[i] = loss
            pbar.set_description(f"Training... Loss: {loss}")
        self.history_ = {}
        self.history_["loss"] = history
        self.module.eval()

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_mudata(
        cls,
        mdata: MuData,
        density_prior_key: Union[
            str, Literal["rna_count_based", "uniform"], None
        ] = "rna_count_based",
        sc_layer: Optional[str] = None,
        sp_layer: Optional[str] = None,
        modalities: Optional[Dict[str, str]] = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
        mdata
            MuData with scRNA and spatial modalities.
        sc_layer
            Layer key in scRNA modality to use for training.
        sp_layer
            Layer key in spatial modality to use for training.
        density_prior_key
            Key in spatial modality obs for density prior.
        modalities
            Mapping from `setup_mudata` param name to modality in mdata.
        """
        setup_method_args = cls._get_setup_method_args(**locals())

        if modalities is None:
            raise ValueError("Modalities cannot be None.")
        modalities = cls._create_modalities_attr_dict(modalities, setup_method_args)

        mudata_fields = [
            fields.MuDataLayerField(
                TANGRAM_REGISTRY_KEYS.SC_KEY,
                sc_layer,
                mod_key=modalities.sc_layer,
                is_count_data=False,
                mod_required=True,
            ),
            fields.MuDataLayerField(
                TANGRAM_REGISTRY_KEYS.SP_KEY,
                sp_layer,
                mod_key=modalities.sp_layer,
                is_count_data=False,
                mod_required=True,
            ),
            fields.MuDataNumericalObsField(
                TANGRAM_REGISTRY_KEYS.DENSITY_KEY,
                density_prior_key,
                mod_key=modalities.density_prior_key,
                required=False,
                mod_required=True,
            ),
        ]
        adata_manager = AnnDataManager(
            fields=mudata_fields,
            setup_method_args=setup_method_args,
            validation_checks=AnnDataManagerValidationCheck(
                check_fully_paired_mudata=False
            ),
        )
        adata_manager.register_fields(mdata, **kwargs)
        sc_state = adata_manager.get_state_registry(TANGRAM_REGISTRY_KEYS.SC_KEY)
        sp_state = adata_manager.get_state_registry(TANGRAM_REGISTRY_KEYS.SP_KEY)
        # Need to access the underlying AnnData field to get these attributes
        if not (
            pd.Index(sc_state[fields.LayerField.COLUMN_NAMES_KEY]).equals(
                sp_state[fields.LayerField.COLUMN_NAMES_KEY]
            ),
        ):
            raise ValueError(
                "The column names of the spatial and single-cell layers must be the same."
            )
        cls.register_manager(adata_manager)

    @classmethod
    def setup_anndata(cls):
        """Not implemented, use `setup_mudata`."""
        raise NotImplementedError(
            "Use `setup_mudata` to setup a MuData object for training."
        )

    def _get_tensor_dict(
        self,
        device: Device,
    ) -> Dict[str, jnp.ndarray]:
        """Get training data for Tangram model.

        Tangram does not minibatch, so we just make a dictionary of
        jnp arrays here.
        """
        tensor_dict = {}
        for key in TANGRAM_REGISTRY_KEYS:
            try:
                tensor_dict[key] = self.adata_manager.get_from_registry(key)
            # When density is missing
            except KeyError:
                continue
            if scipy.sparse.issparse(tensor_dict[key]):
                tensor_dict[key] = tensor_dict[key].toarray()
            elif isinstance(tensor_dict[key], pd.DataFrame):
                tensor_dict[key] = tensor_dict[key].values
            else:
                tensor_dict[key] = tensor_dict[key]
            tensor_dict[key] = _asarray(tensor_dict[key], device=device)

        return tensor_dict

    @staticmethod
    def project_cell_annotations(
        adata_sc: AnnData, adata_sp: AnnData, mapper: np.ndarray, labels: pd.Series
    ) -> pd.DataFrame:
        """Project cell annotations to spatial data.

        Parameters
        ----------
        adata_sc
            AnnData object with single-cell data.
        adata_sp
            AnnData object with spatial data.
        mapper
            Mapping from single-cell to spatial data.
        labels
            Cell annotations to project.

        Returns
        -------
        Projected annotations as a :class:`pd.DataFrame` with shape (n_sp, n_labels).
        """
        if len(labels) != adata_sc.shape[0]:
            raise ValueError(
                "The number of labels must match the number of cells in the sc AnnData object."
            )
        cell_type_df = pd.get_dummies(labels)
        projection = mapper.T @ cell_type_df.values
        return pd.DataFrame(
            index=adata_sp.obs_names, columns=cell_type_df.columns, data=projection
        )

    @staticmethod
    def project_genes(
        adata_sc: AnnData, adata_sp: AnnData, mapper: np.ndarray
    ) -> AnnData:
        """Project gene expression to spatial data.

        Parameters
        ----------
        adata_sc
            AnnData object with single-cell data.
        adata_sp
            AnnData object with spatial data.
        mapper
            Mapping from single-cell to spatial data.

        Returns
        -------
        :class:`anndata.AnnData` object with projected gene expression.
        """
        adata_ge = AnnData(
            X=mapper.T @ adata_sc.X,
            obs=adata_sp.obs,
            var=adata_sc.var,
            uns=adata_sc.uns,
        )
        return adata_ge
