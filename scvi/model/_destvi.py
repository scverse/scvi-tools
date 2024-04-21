from __future__ import annotations

import logging
from collections import OrderedDict
from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import torch
from anndata import AnnData

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import (
    CategoricalJointObsField,
    LayerField,
    NumericalObsField,
    CategoricalObsField,
)
from scvi.model import CondSCVI
from scvi.model.base import BaseModelClass, UnsupervisedTrainingMixin
from scvi.module import MRDeconv
from scvi.utils import setup_anndata_dsp
from scvi.utils._docstrings import devices_dsp

logger = logging.getLogger(__name__)


class DestVI(UnsupervisedTrainingMixin, BaseModelClass):
    """Multi-resolution deconvolution of Spatial Transcriptomics data (DestVI) :cite:p:`Lopez22`. Most users will use the alternate constructor (see example).

    Parameters
    ----------
    st_adata
        spatial transcriptomics AnnData object that has been registered via :meth:`~scvi.model.DestVI.setup_anndata`.
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
        Keyword args for :class:`~scvi.modules.MRDeconv`

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
    2. :doc:`/tutorials/notebooks/spatial/DestVI_in_R`
    """

    _module_cls = MRDeconv

    def __init__(
        self,
        st_adata: AnnData,
        cell_type_mapping: np.ndarray,
        decoder_state_dict: OrderedDict,
        px_decoder_state_dict: OrderedDict,
        px_r: np.ndarray,
        n_hidden: int,
        n_latent: int,
        n_layers: int,
        n_batch_sc: int,
        dropout_decoder: float,
        l1_reg: float,
        sc_batch_mapping: list[str],
        **module_kwargs,
    ):
        super().__init__(st_adata)
        if sc_batch_mapping is not None:
            st_sc_batch_mapping = self.adata_manager.get_state_registry('batch_index_sc')['categorical_mapping']
            assert set(st_sc_batch_mapping).issubset(set(sc_batch_mapping)), (               
                f'Spatial model has other covariates than single cell model, {set(st_sc_batch_mapping) - set(sc_batch_mapping)}'
            )

        self.module = self._module_cls(
            n_spots=st_adata.n_obs,
            n_labels=cell_type_mapping.shape[0],
            n_batch=self.summary_stats.n_batch,
            n_batch_sc=n_batch_sc,
            decoder_state_dict=decoder_state_dict,
            px_decoder_state_dict=px_decoder_state_dict,
            px_r=px_r,
            n_genes=st_adata.n_vars,
            n_latent=n_latent,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_decoder=dropout_decoder,
            l1_reg=l1_reg,
            **module_kwargs,
        )
        self.cell_type_mapping = cell_type_mapping
        self._model_summary_string = "DestVI Model"
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    def from_rna_model(
        cls,
        st_adata: AnnData,
        sc_model: CondSCVI,
        vamp_prior_p: int = 15,
        l1_reg: float = 0.0,
        **module_kwargs,
    ):
        """Alternate constructor for exploiting a pre-trained model on a RNA-seq dataset.

        Parameters
        ----------
        st_adata
            registered anndata object
        sc_model
            trained CondSCVI model
        vamp_prior_p
            number of mixture parameter for VampPrior calculations
        l1_reg
            Scalar parameter indicating the strength of L1 regularization on cell type proportions.
            A value of 50 leads to sparser results.
        **model_kwargs
            Keyword args for :class:`~scvi.model.DestVI`
        """
        decoder_state_dict = sc_model.module.decoder.state_dict()
        px_decoder_state_dict = sc_model.module.px_decoder.state_dict()
        px_r = sc_model.module.px_r.detach().cpu().numpy()
        mapping = sc_model.adata_manager.get_state_registry(
            REGISTRY_KEYS.LABELS_KEY
        ).categorical_mapping

        dropout_decoder = sc_model.module.dropout_rate
        if vamp_prior_p is None:
            mean_vprior = None
            var_vprior = None
        else:
            mean_vprior, var_vprior, mp_vprior = sc_model.get_vamp_prior(
                sc_model.adata, p=vamp_prior_p
            ).values()

        sc_batch_mapping = (
            sc_model.adata_manager.get_state_registry(
                REGISTRY_KEYS.BATCH_KEY
            )
        )['categorical_mapping']

        return cls(
            st_adata,
            mapping,
            decoder_state_dict,
            px_decoder_state_dict,
            px_r,
            sc_model.module.n_hidden,
            sc_model.module.n_latent,
            sc_model.module.n_layers,
            sc_model.module.n_batch,
            mean_vprior=mean_vprior,
            var_vprior=var_vprior,
            mp_vprior=mp_vprior,
            dropout_decoder=dropout_decoder,
            l1_reg=l1_reg,
            sc_batch_mapping=sc_batch_mapping,
            **module_kwargs,
        )

    @torch.inference_mode()
    def get_proportions(
        self,
        keep_noise: bool = False,
        normalize: bool = True,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
    ) -> pd.DataFrame:
        """Returns the estimated cell type proportion for the spatial data.

        Shape is n_cells x n_labels OR n_cells x (n_labels + 1) if keep_noise.

        Parameters
        ----------
        keep_noise
            whether to account for the noise term as a standalone cell type in the proportion estimate.
        normalize
            whether to normalize the proportions to sum to 1.
        indices
            Indices of cells in adata to use. Only used if amortization. If `None`, all cells are used.
        batch_size
            Minibatch size for data loading into model. Only used if amortization. Defaults to `scvi.settings.batch_size`.
        """
        self._check_if_trained()

        column_names = self.cell_type_mapping
        index_names = self.adata.obs.index
        if keep_noise:
            column_names = np.append(column_names, "noise_term")

        if self.module.amortization in ["both", "proportion"]:
            stdl = self._make_data_loader(
                adata=self.adata, indices=indices, batch_size=batch_size
            )
            prop_ = []
            for tensors in stdl:
                inference_inputs = self.module._get_inference_input(tensors)
                outputs = self.module.inference(**inference_inputs)
                generative_inputs = self.module._get_generative_input(tensors, outputs)
                prop_local = self.module.generative(**generative_inputs)["v"]
                prop_ += [prop_local.cpu()]
            data = torch.cat(prop_).numpy()
            if indices:
                index_names = index_names[indices]
        else:
            if indices is not None:
                logger.info(
                    "No amortization for proportions, ignoring indices and returning results for the full data"
                )
            data = torch.nn.functional.softplus(self.module.V).transpose(1, 0).detach().cpu().numpy()
        if normalize:
            data = data / data.sum(axis=1, keepdims=True)
        if not keep_noise:
            data = data[:, :-1]

        return pd.DataFrame(
            data=data,
            columns=column_names,
            index=index_names,
        )
        
    @torch.inference_mode()
    def get_fine_celltypes(
        self,
        sc_model: CondSCVI,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        return_numpy: bool = False,
    ) -> np.ndarray | dict[str, pd.DataFrame]:
        """Returns the estimated cell-type specific latent space for the spatial data.

        Parameters
        ----------
        sc_model
            trained CondSCVI model
        indices
            Indices of cells in adata to use. Only used if amortization. If `None`, all cells are used.
        batch_size
            Minibatch size for data loading into model. Only used if amortization. Defaults to `scvi.settings.batch_size`.
        return_numpy
            if activated, will return a numpy array of shape is n_spots x n_latent x n_labels.
        """
        self._check_if_trained()

        column_names = [str(i) for i in np.arange(self.module.n_latent)]
        index_names = self.adata.obs.index

        if self.module.amortization in ["both", "latent"]:
            stdl = self._make_data_loader(
                adata=self.adata, indices=indices, batch_size=batch_size
            )
            gamma_ = []
            proportions_modes_ = []
            for tensors in stdl:
                inference_inputs = self.module._get_inference_input(tensors)
                outputs = self.module.inference(**inference_inputs)
                generative_inputs = self.module._get_generative_input(tensors, outputs)
                generative_outputs = self.module.generative(**generative_inputs)
                gamma_local = generative_outputs["gamma"]
                if self.module.prior_mode == 'mog':
                    proportions_modes_local = generative_outputs['proportion_modes'] # pmc
                    gamma_local = gamma_local # pncm
                else:
                    proportions_modes_local = torch.ones(gamma_local.shape[0], 1, 1)
                    gamma_local = gamma_local.squeeze(0) # pncm
                gamma_ += [gamma_local.cpu()]
                proportions_modes_ += [proportions_modes_local.cpu()]
                
            proportions_modes = torch.cat(proportions_modes_, dim=-1).numpy()
            gamma = torch.cat(gamma_, dim=-1).numpy()
        else:
            if indices is not None:
                logger.info(
                    "No amortization for latent values, ignoring adata and returning results for the full data"
                )
            gamma = self.module.gamma.detach().cpu().numpy()
        
        sc_latent_distribution = sc_model.get_latent_representation(return_dist=True)


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
            Indices of cells in adata to use. Only used if amortization. If `None`, all cells are used.
        batch_size
            Minibatch size for data loading into model. Only used if amortization. Defaults to `scvi.settings.batch_size`.
        return_numpy
            if activated, will return a numpy array of shape is n_spots x n_latent x n_labels.
        """
        self._check_if_trained()

        column_names = [str(i) for i in np.arange(self.module.n_latent)]
        index_names = self.adata.obs.index

        if self.module.amortization in ["both", "latent"]:
            stdl = self._make_data_loader(
                adata=self.adata, indices=indices, batch_size=batch_size
            )
            gamma_ = []
            for tensors in stdl:
                inference_inputs = self.module._get_inference_input(tensors)
                outputs = self.module.inference(**inference_inputs)
                generative_inputs = self.module._get_generative_input(tensors, outputs)
                generative_outputs = self.module.generative(**generative_inputs)
                gamma_local = generative_outputs["gamma"]
                if self.module.prior_mode == 'mog':
                    proportions_model_local = generative_outputs['proportion_modes']
                    gamma_local = torch.einsum('pncm,pmc->ncm', gamma_local, proportions_model_local)
                else:
                    gamma_local = gamma_local.squeeze(0)
                gamma_ += [gamma_local.cpu()]
            data = torch.cat(gamma_, dim=-1).numpy()
            if indices is not None:
                index_names = index_names[indices]
        else:
            if indices is not None:
                logger.info(
                    "No amortization for latent values, ignoring adata and returning results for the full data"
                )
            data = self.module.gamma.detach().cpu().numpy()

        data = np.transpose(data, (2, 0, 1))
        if return_numpy:
            return data
        else:
            res = {}
            for i, ct in enumerate(self.cell_type_mapping):
                res[ct] = pd.DataFrame(
                    data=data[:, :, i], columns=column_names, index=index_names
                )
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
    ) -> np.ndarray (np.ndarray, np.ndarray):
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
            For distributions with no closed-form mean (e.g., `logistic normal`), how many Monte Carlo
            samples to take for computing mean.
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
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )
        latent = []
        latent_qzm = []
        latent_qzv = []
        for tensors in scdl:
            inference_inputs = self.module._get_inference_input(tensors)
            z, qz, _ = self.module.inference(**inference_inputs, n_samples=mc_samples).values()
            if give_mean:
                latent += [qz.loc.cpu()]
            else:
                latent += [z.cpu()]
            latent_qzm += [qz.loc.cpu()]
            latent_qzv += [qz.scale.square().cpu()]
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
        
        cell_type_mapping_extended = list(self.cell_type_mapping) + ['noise']

        if label not in cell_type_mapping_extended:
            raise ValueError("Unknown cell type")
        y = cell_type_mapping_extended.index(label)

        stdl = self._make_data_loader(self.adata, indices=indices, batch_size=batch_size)
        scale = []
        for tensors in stdl:
            inference_inputs = self.module._get_inference_input(tensors)
            outputs = self.module.inference(**inference_inputs)
            generative_inputs = self.module._get_generative_input(tensors, outputs)
            px_scale = self.module.generative(**generative_inputs)["px_mu"][:, y, :]

            scale += [px_scale.cpu()]

        data = torch.cat(scale).numpy()
        column_names = self.adata.var.index
        index_names = self.adata.obs.index
        if indices is not None:
            index_names = index_names[indices]
        return pd.DataFrame(data=data, columns=column_names, index=index_names)
    
    @torch.inference_mode()
    def get_expression_for_ct(
        self,
        label: str,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        return_sparse_array: bool = False,
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
        return_sparse_array
            If `True`, returns a sparse array instead of a dataframe.

        Returns
        -------
        Pandas dataframe of gene_expression
        """
        self._check_if_trained()
        cell_type_mapping_extended = list(self.cell_type_mapping) + ['noise']

        if label not in cell_type_mapping_extended:
            raise ValueError("Unknown cell type")
        y = cell_type_mapping_extended.index(label)

        stdl = self._make_data_loader(self.adata, indices=indices, batch_size=batch_size)
        expression_ct = []
        for tensors in stdl:
            inference_inputs = self.module._get_inference_input(tensors)
            outputs = self.module.inference(**inference_inputs)
            generative_inputs = self.module._get_generative_input(tensors, outputs)
            generative_outputs = self.module.generative(**generative_inputs)
            px_scale, proportions = generative_outputs['px_mu'], generative_outputs['v']
            px_scale = torch.einsum('mkl,mk->mkl', px_scale, proportions)
            px_scale_proportions = px_scale[:, y, :]/px_scale.sum(dim=1)
            x_ct = inference_inputs['x'].to(px_scale_proportions.device) * px_scale_proportions
            expression_ct += [x_ct.cpu()]

        data = torch.cat(expression_ct).numpy()
        if return_sparse_array:
            data = csr_matrix(data.T)
            return data
        else:
            column_names = self.adata.var.index
            index_names = self.adata.obs.index
            if indices is not None:
                index_names = index_names[indices]
            return pd.DataFrame(data=data, columns=column_names, index=index_names)

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int = 2000,
        lr: float = 0.003,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float = 1.0,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        n_epochs_kl_warmup: int = 200,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **kwargs,
    ):
        """Trains the model using MAP inference.

        Parameters
        ----------
        max_epochs
            Number of epochs to train for
        lr
            Learning rate for optimization.
        %(param_accelerator)s
        %(param_devices)s
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        shuffle_set_split
            Whether to shuffle indices before splitting. If `False`, the val, train, and test set are split in the
            sequential order of the data according to `validation_size` and `train_size` percentages.
        batch_size
            Minibatch size to use during training.
        n_epochs_kl_warmup
            number of epochs needed to reach unit kl weight in the elbo
        datasplitter_kwargs
            Additional keyword arguments passed into :class:`~scvi.dataloaders.DataSplitter`.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        update_dict = {
            "lr": lr,
            "n_epochs_kl_warmup": n_epochs_kl_warmup,
        }
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict
        super().train(
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            batch_size=batch_size,
            datasplitter_kwargs=datasplitter_kwargs,
            plan_kwargs=plan_kwargs,
            **kwargs,
        )

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        sc_batch_key: str | None = None,
        categorical_covariate_keys: Sequence[str] | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        %(param_batch_key)s
        sc_batch_key:
        Categorical covariate keys need to line up with single cell model.
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        # add index for each cell (provided to pyro plate for correct minibatching)
        adata.obs["_indices"] = np.arange(adata.n_obs)
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            NumericalObsField(REGISTRY_KEYS.INDICES_KEY, "_indices"),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField("batch_index_sc", sc_batch_key),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
