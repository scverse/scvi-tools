from __future__ import annotations

import logging
from collections import OrderedDict
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import torch

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data._constants import _SETUP_ARGS_KEY
from scvi.data.fields import (
    CategoricalObsField,
    LayerField,
    NumericalObsField,
)
from scvi.model.base import BaseModelClass, UnsupervisedTrainingMixin
from scvi.model.base._archesmixin import _get_loaded_data
from scvi.module import MRDeconv
from scvi.utils import setup_anndata_dsp

if TYPE_CHECKING:
    from collections.abc import Sequence

    from anndata import AnnData

    from scvi.model import CondSCVI

logger = logging.getLogger(__name__)


class DestVI(UnsupervisedTrainingMixin, BaseModelClass):
    """Multi-resolution deconvolution of Spatial Transcriptomics data (DestVI) :cite:p:`Lopez22`.

    Most users will use the alternate constructor (see example).

    Parameters
    ----------
    st_adata
        spatial transcriptomics AnnData object that has been registered via
        :meth:`~scvi.model.DestVI.setup_anndata`.
    cell_type_mapping
        mapping between numerals and cell type labels
    decoder_state_dict
        state_dict from the decoder of the CondSCVI model
    px_decoder_state_dict
        state_dict from the px_decoder of the CondSCVI model
    px_r
        parameters for the px_r tensor in the CondSCVI model
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    **module_kwargs
        Keyword args for :class:`~scvi.module.MRDeconv`

    Examples
    --------
    >>> sc_adata = anndata.read_h5ad(path_to_scRNA_anndata)
    >>> scvi.model.CondSCVI.setup_anndata(sc_adata)
    >>> sc_model = scvi.model.CondSCVI(sc_adata)
    >>> st_adata = anndata.read_h5ad(path_to_ST_anndata)
    >>> DestVI.setup_anndata(st_adata)
    >>> spatial_model = DestVI.from_rna_model(st_adata, sc_model)
    >>> spatial_model.train(max_epochs=2000)
    >>> st_adata.obsm["proportions"] = spatial_model.get_proportions(st_adata)
    >>> gamma = spatial_model.get_gamma(st_adata)

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/spatial/DestVI_tutorial`
    2. :doc:`/tutorials/notebooks/r/DestVI_in_R`
    """

    _module_cls = MRDeconv

    def __init__(
        self,
        st_adata: AnnData,
        cell_type_mapping: np.ndarray,
        decoder_state_dict: OrderedDict,
        px_decoder_state_dict: OrderedDict,
        px_r: torch.tensor,
        per_ct_bias: torch.tensor,
        n_hidden: int,
        n_latent: int,
        n_layers: int,
        dropout_decoder: float,
        **module_kwargs,
    ):
        super().__init__(st_adata)

        self.module = self._module_cls(
            n_spots=st_adata.n_obs,
            n_labels=cell_type_mapping.shape[0],
            n_batch=self.summary_stats.n_batch,
            decoder_state_dict=decoder_state_dict,
            px_decoder_state_dict=px_decoder_state_dict,
            px_r=px_r,
            per_ct_bias=per_ct_bias,
            n_genes=st_adata.n_vars,
            n_latent=n_latent,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_decoder=dropout_decoder,
            **module_kwargs,
        )
        self.cell_type_mapping = cell_type_mapping
        self.cell_type_mapping_extended = list(self.cell_type_mapping) + [
            f"additional_{i}" for i in range(self.module.add_celltypes)
        ]
        self._model_summary_string = "DestVI Model"
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    def from_rna_model(
        cls,
        st_adata: AnnData,
        sc_model: CondSCVI,
        vamp_prior_p: int = 15,
        anndata_setup_kwargs: dict | None = None,
        **module_kwargs,
    ):
        """Alternate constructor for exploiting a pre-trained model on a RNA-seq dataset.

        Parameters
        ----------
        st_adata
            anndata object will be registered
        sc_model
            trained CondSCVI model or path to a trained model
        vamp_prior_p
            number of mixture parameter for VampPrior calculations
        anndata_setup_kwargs
            Keyword args for :meth:`~scvi.model.DestVI.setup_anndata`
        **module_kwargs
            Keyword args for :class:`~scvi.model.MRDeconv`
        """
        attr_dict, var_names, load_state_dict, _ = _get_loaded_data(sc_model)
        registry = attr_dict.pop("registry_")

        decoder_state_dict = OrderedDict(
            (i[8:], load_state_dict[i])
            for i in load_state_dict.keys()
            if i.split(".")[0] == "decoder"
        )
        px_decoder_state_dict = OrderedDict(
            (i[11:], load_state_dict[i])
            for i in load_state_dict.keys()
            if i.split(".")[0] == "px_decoder"
        )
        px_r = load_state_dict["px_r"]
        per_ct_bias = load_state_dict["per_ct_bias"]
        mapping = registry["field_registries"]["labels"]["state_registry"]["categorical_mapping"]

        dropout_decoder = attr_dict["init_params_"]["non_kwargs"]["dropout_rate"]
        if vamp_prior_p is None:
            mean_vprior = None
            var_vprior = None
        elif attr_dict["init_params_"]["kwargs"]["module_kwargs"]["prior"] == "mog":
            mean_vprior = load_state_dict["prior_means"].clone().detach()
            var_vprior = torch.exp(load_state_dict["prior_log_std"]) ** 2
            mp_vprior = torch.nn.Softmax(dim=-1)(load_state_dict["prior_logits"])
        else:
            assert sc_model is not str, (
                "VampPrior requires loading CondSCVI model and providing it"
            )
            mean_vprior, var_vprior, mp_vprior = sc_model.get_vamp_prior(
                sc_model.adata, p=vamp_prior_p
            ).values()

        if anndata_setup_kwargs is None:
            anndata_setup_kwargs = {}

        cls.setup_anndata(
            st_adata,
            source_registry=registry,
            extend_categories=True,
            **anndata_setup_kwargs,
            **registry[_SETUP_ARGS_KEY],
        )

        return cls(
            st_adata,
            mapping,
            decoder_state_dict,
            px_decoder_state_dict,
            px_r,
            per_ct_bias,
            sc_model.module.n_hidden,
            sc_model.module.n_latent,
            sc_model.module.n_layers,
            mean_vprior=mean_vprior,
            var_vprior=var_vprior,
            mp_vprior=mp_vprior,
            dropout_decoder=dropout_decoder,
            **module_kwargs,
        )

    @torch.inference_mode()
    def get_proportions(
        self,
        keep_additional: bool = False,
        normalize: bool = True,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
    ) -> pd.DataFrame:
        """Returns the estimated cell type proportion for the spatial data.

        Shape is n_cells x n_labels OR n_cells x (n_labels + add_celltypes) if keep_additional.

        Parameters
        ----------
        keep_additional
            whether to account for the additional cell-types as standalone cell types
            in the proportion estimate.
        normalize
            whether to normalize the proportions to sum to 1.
        indices
            Indices of cells in adata to use. Only used if amortization. If `None`, all cells are
            used.
        batch_size
            Minibatch size for data loading into model. Only used if amortization. Defaults to
            `scvi.settings.batch_size`.
        """
        self._check_if_trained()

        column_names = self.cell_type_mapping
        index_names = self.adata.obs.index
        if keep_additional:
            column_names = list(self.cell_type_mapping_extended)
        else:
            column_names = list(self.cell_type_mapping)

        if self.module.amortization in ["both", "proportion"]:
            stdl = self._make_data_loader(adata=self.adata, indices=indices, batch_size=batch_size)
            prop_ = []
            for tensors in stdl:
                inference_inputs = self.module._get_inference_input(tensors)
                outputs = self.module.inference(**inference_inputs)
                generative_inputs = self.module._get_generative_input(tensors, outputs)
                prop_local = self.module.generative(**generative_inputs)["v"][0, ...]
                prop_ += [prop_local.cpu()]
            data = torch.cat(prop_).detach().numpy()
            if indices:
                index_names = index_names[indices]
        else:
            data = (
                torch.nn.functional.softplus(self.module.V).transpose(0, 1).detach().cpu().numpy()
            )
        if not keep_additional:
            data = data[:, : -self.module.add_celltypes]
        if normalize:
            data = data / data.sum(axis=1, keepdims=True)

        return pd.DataFrame(
            data=data,
            columns=column_names,
            index=index_names,
        )

    @torch.inference_mode()
    def get_gamma(
        self,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        return_numpy: bool = False,
    ) -> np.ndarray | dict[str, pd.DataFrame]:
        """Returns the estimated cell-type specific latent space for the spatial data.

        Parameters
        ----------
        indices
            Indices of cells in adata to use. Only used if amortization. If `None`, all cells are
            used.
        batch_size
            Minibatch size for data loading into model. Only used if amortization. Defaults to
            `scvi.settings.batch_size`.
        return_numpy
            if activated, will return a numpy array of shape is n_spots x n_latent x n_labels.
        """
        self._check_if_trained()

        column_names = [str(i) for i in np.arange(self.module.n_latent)]
        index_names = self.adata.obs.index

        if self.module.amortization in ["both", "latent"]:
            stdl = self._make_data_loader(adata=self.adata, indices=indices, batch_size=batch_size)
            gamma_ = []
            for tensors in stdl:
                inference_inputs = self.module._get_inference_input(tensors)
                outputs = self.module.inference(**inference_inputs)
                generative_inputs = self.module._get_generative_input(tensors, outputs)
                generative_outputs = self.module.generative(**generative_inputs)
                gamma_local = generative_outputs["gamma"][0, ...]
                if self.module.prior_mode == "mog":
                    proportions_model_local = generative_outputs["proportion_modes"][0, ...]
                    gamma_local = torch.einsum(
                        "pncm,pmc->ncm", gamma_local, proportions_model_local
                    )
                else:
                    gamma_local = gamma_local[0, ...].squeeze(0)
                gamma_ += [gamma_local.cpu()]
            data = torch.cat(gamma_, dim=-1).numpy()
            if indices is not None:
                index_names = index_names[indices]
        else:
            data = self.module.gamma.detach().cpu().numpy()
        data = np.transpose(data, (2, 0, 1))
        if return_numpy:
            return data
        else:
            res = {}
            for i, ct in enumerate(self.cell_type_mapping):
                res[ct] = pd.DataFrame(data=data[:, :, i], columns=column_names, index=index_names)
            return res

    @torch.inference_mode()
    def get_latent_representation(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        give_mean: bool = True,
        mc_samples: int = 5000,
        batch_size: int | None = None,
        return_dist: bool = False,
    ) -> np.ndarray(np.ndarray, np.ndarray):
        """Return the latent representation for each cell.

        This is typically denoted as :math:`z_n`.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Give mean of distribution or sample from it.
        mc_samples
            For distributions with no closed-form mean (e.g., `logistic normal`),
            how many Monte Carlo samples to take for computing mean.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_dist
            Return (mean, variance) of distributions instead of just the mean.
            If `True`, ignores `give_mean` and `mc_samples`. In the case of the latter,
            `mc_samples` is used to compute the mean of a transformed distribution.
            If `return_dist` is true the untransformed mean and variance are returned.

        Returns
        -------
        Low-dimensional representation for each cell or a tuple containing its mean and variance.
        """
        assert self.module.n_latent_amortization is not None, (
            "Model has no latent representation for amortized values."
        )
        self._check_if_trained(warn=False)

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        latent = []
        latent_qzm = []
        latent_qzv = []
        for tensors in scdl:
            inference_inputs = self.module._get_inference_input(tensors)
            inference_outputs = self.module.inference(**inference_inputs, n_samples=mc_samples)
            z = inference_outputs["z"][0, ...]
            qz = inference_outputs["qz"]
            if give_mean:
                latent += [qz.loc[0, ...].cpu()]
            else:
                latent += [z.cpu()]
            latent_qzm += [qz.loc[0, ...].cpu()]
            latent_qzv += [qz.scale[0, ...].square().cpu()]
        return (
            (torch.cat(latent_qzm).numpy(), torch.cat(latent_qzv).numpy())
            if return_dist
            else torch.cat(latent).numpy()
        )

    @torch.inference_mode()
    def get_scale_for_ct(
        self,
        label: str,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
    ) -> pd.DataFrame:
        r"""Return the scaled parameter of the NB for every spot in queried cell types.

        Parameters
        ----------
        label
            cell type of interest
        indices
            Indices of cells in self.adata to use. If `None`, all cells are used.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        Pandas dataframe of gene_expression
        """
        self._check_if_trained()
        self._validate_anndata()

        cell_type_mapping_extended = list(self.cell_type_mapping) + [
            f"additional_{i}" for i in range(self.module.add_celltypes)
        ]

        if label not in cell_type_mapping_extended:
            raise ValueError("Unknown cell type")
        y = cell_type_mapping_extended.index(label)

        stdl = self._make_data_loader(self.adata, indices=indices, batch_size=batch_size)
        scale = []
        for tensors in stdl:
            inference_inputs = self.module._get_inference_input(tensors)
            outputs = self.module.inference(**inference_inputs)
            generative_inputs = self.module._get_generative_input(tensors, outputs)
            px_scale = self.module.generative(**generative_inputs)["px_mu"][0, :, y, :]

            scale += [px_scale.cpu()]

        data = torch.cat(scale).numpy()
        column_names = self.adata.var.index
        index_names = self.adata.obs.index
        if indices is not None:
            index_names = index_names[indices]
        return pd.DataFrame(data=data, columns=column_names, index=index_names)

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        smoothed_layer: str | None = None,
        batch_key: str | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        smoothed_layer
            param that...
        %(param_batch_key)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        # add index for each cell (provided to pyro plate for correct minibatching)
        adata.obs["_indices"] = np.arange(adata.n_obs)
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            NumericalObsField(REGISTRY_KEYS.INDICES_KEY, "_indices"),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
        ]
        if smoothed_layer is not None:
            anndata_fields.append(LayerField("x_smoothed", smoothed_layer, is_count_data=True))
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
