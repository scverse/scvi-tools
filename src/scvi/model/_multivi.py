from __future__ import annotations

import logging
import warnings
from collections.abc import Iterable as IterableClass
from functools import partial
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import torch
from mudata import MuData
from scipy.sparse import csr_matrix, vstack
from torch.distributions import Normal

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager, fields
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
    ProteinObsmField,
)
from scvi.model._utils import (
    _get_batch_code_from_category,
    scatac_raw_counts_properties,
    scrna_raw_counts_properties,
    use_distributed_sampler,
)
from scvi.model.base import (
    ArchesMixin,
    BaseModelClass,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.model.base._de_core import _de_core
from scvi.module import MULTIVAE
from scvi.train import AdversarialTrainingPlan
from scvi.utils import track
from scvi.utils._docstrings import de_dsp, devices_dsp, setup_anndata_dsp

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from typing import Literal

    from anndata import AnnData

    from scvi._types import AnnOrMuData, Number

logger = logging.getLogger(__name__)


class MULTIVI(VAEMixin, UnsupervisedTrainingMixin, BaseModelClass, ArchesMixin):
    """Integration of multi-modal and single-modality data :cite:p:`AshuachGabitto21`.

    MultiVI is used to integrate multiomic datasets with single-modality (expression
    or accessibility) datasets.

    Parameters
    ----------
    adata
        AnnData/MuData object that has been registered via
        :meth:`~scvi.model.MULTIVI.setup_anndata` or :meth:`~scvi.model.MULTIVI.setup_mudata`.
    n_genes
        The number of gene expression features (genes).
    n_regions
        The number of accessibility features (genomic regions).
    modality_weights
        Weighting scheme across modalities. One of the following:
        * ``"equal"``: Equal weight in each modality
        * ``"universal"``: Learn weights across modalities w_m.
        * ``"cell"``: Learn weights across modalities and cells. w_{m,c}
    modality_penalty
        Training Penalty across modalities. One of the following:
        * ``"Jeffreys"``: Jeffreys penalty to align modalities
        * ``"MMD"``: MMD penalty to align modalities
        * ``"None"``: No penalty
    n_hidden
        Number of nodes per hidden layer. If `None`, defaults to square root
        of number of regions.
    n_latent
        Dimensionality of the latent space. If `None`, defaults to square root
        of `n_hidden`.
    n_layers_encoder
        Number of hidden layers used for encoder NNs.
    n_layers_decoder
        Number of hidden layers used for decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    model_depth
        Model sequencing depth / library size.
    region_factors
        Include region-specific factors in the model.
    gene_dispersion
        One of the following
        * ``'gene'`` - genes_dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - genes_dispersion can differ between different batches
        * ``'gene-label'`` - genes_dispersion can differ between different labels
    protein_dispersion
        One of the following
        * ``'protein'`` - protein_dispersion parameter is constant per protein across cells
        * ``'protein-batch'`` - protein_dispersion can differ between different batches NOT TESTED
        * ``'protein-label'`` - protein_dispersion can differ between different labels NOT TESTED
    latent_distribution
        One of
        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
    deeply_inject_covariates
        Whether to deeply inject covariates into all layers of the decoder. If False,
        covariates will only be included in the input layer.
    fully_paired
        allows the simplification of the model if the data is fully paired. Currently ignored.
    **model_kwargs
        Keyword args for :class:`~scvi.module.MULTIVAE`

    Examples
    --------
    >>> adata_rna = anndata.read_h5ad(path_to_rna_anndata)
    >>> adata_atac = scvi.data.read_10x_atac(path_to_atac_anndata)
    >>> adata_protein = anndata.read_h5ad(path_to_protein_anndata)
    >>> mdata = MuData({"rna": adata_rna, "protein": adata_protein, "atac": adata_atac})
    >>> scvi.model.MULTIVI.setup_mudata(mdata, batch_key="batch",
    >>> modalities={"rna_layer": "rna", "protein_layer": "protein", "batch_key": "rna",
    >>>             "atac_layer": "atac"})
    >>> vae = scvi.model.MULTIVI(mdata)
    >>> vae.train()

    Notes (for using setup_anndata)
    -----
    * The model assumes that the features are organized so that all expression features are
       consecutive, followed by all accessibility features. For example, if the data has 100 genes
       and 250 genomic regions, the model assumes that the first 100 features are genes, and the
       next 250 are the regions.

    * The main batch annotation, specified in ``setup_anndata``, should correspond to
       the modality each cell originated from. This allows the model to focus mixing efforts, using
       an adversarial component, on mixing the modalities. Other covariates can be specified using
       the `categorical_covariate_keys` argument.
    """

    _module_cls = MULTIVAE
    _training_plan_cls = AdversarialTrainingPlan

    def __init__(
        self,
        adata: AnnOrMuData,
        n_genes: int | None = None,
        n_regions: int | None = None,
        modality_weights: Literal["equal", "cell", "universal"] = "equal",
        modality_penalty: Literal["Jeffreys", "MMD", "None"] = "Jeffreys",
        n_hidden: int | None = None,
        n_latent: int | None = None,
        n_layers_encoder: int = 2,
        n_layers_decoder: int = 2,
        dropout_rate: float = 0.1,
        region_factors: bool = True,
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb",
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        latent_distribution: Literal["normal", "ln"] = "normal",
        deeply_inject_covariates: bool = False,
        encode_covariates: bool = False,
        fully_paired: bool = False,
        protein_dispersion: Literal["protein", "protein-batch", "protein-label"] = "protein",
        **model_kwargs,
    ):
        super().__init__(adata)
        if "n_proteins" in self.summary_stats:
            self.protein_state_registry = self.adata_manager.get_state_registry(
                REGISTRY_KEYS.PROTEIN_EXP_KEY
            )

        if n_genes is None or n_regions is None:
            assert isinstance(adata, MuData), (
                "n_genes and n_regions must be provided if using AnnData"
            )
            n_genes = self.summary_stats.get("n_vars", 0)
            n_regions = self.summary_stats.get("n_atac", 0)

        prior_mean, prior_scale = None, None
        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else []
        )

        use_size_factor_key = REGISTRY_KEYS.SIZE_FACTOR_KEY in self.adata_manager.data_registry

        if "n_proteins" in self.summary_stats:
            n_proteins = self.summary_stats.n_proteins
        else:
            n_proteins = 0

        self.module = self._module_cls(
            n_input_genes=n_genes,
            n_input_regions=n_regions,
            n_input_proteins=n_proteins,
            modality_weights=modality_weights,
            modality_penalty=modality_penalty,
            n_batch=self.summary_stats.n_batch,
            n_obs=adata.n_obs,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers_encoder=n_layers_encoder,
            n_layers_decoder=n_layers_decoder,
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            n_cats_per_cov=n_cats_per_cov,
            dropout_rate=dropout_rate,
            region_factors=region_factors,
            gene_likelihood=gene_likelihood,
            gene_dispersion=dispersion,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            use_size_factor_key=use_size_factor_key,
            latent_distribution=latent_distribution,
            deeply_inject_covariates=deeply_inject_covariates,
            encode_covariates=encode_covariates,
            protein_background_prior_mean=prior_mean,
            protein_background_prior_scale=prior_scale,
            protein_dispersion=protein_dispersion,
            **model_kwargs,
        )
        self._model_summary_string = (
            f"MultiVI Model with the following params: \nn_genes: {n_genes}, "
            f"n_regions: {n_regions}, n_proteins: {n_proteins}, n_hidden: {self.module.n_hidden}, "
            f"n_latent: {self.module.n_latent}, n_layers_encoder: {n_layers_encoder}, "
            f"n_layers_decoder: {n_layers_decoder}, dropout_rate: {dropout_rate}, "
            f"latent_distribution: {latent_distribution}, "
            f"deep injection: {deeply_inject_covariates}, gene_likelihood: {gene_likelihood}, "
            f"gene_dispersion:{dispersion}, Mod.Weights: {modality_weights}, "
            f"Mod.Penalty: {modality_penalty}, protein_dispersion: {protein_dispersion}"
        )
        self.fully_paired = fully_paired
        self.n_latent = n_latent
        self.init_params_ = self._get_init_params(locals())
        self.n_genes = n_genes
        self.n_regions = n_regions
        self.n_proteins = n_proteins
        self.get_normalized_function_name = "get_normalized_accessibility"

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int = 500,
        lr: float = 1e-4,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float | None = None,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        weight_decay: float = 1e-3,
        eps: float = 1e-08,
        early_stopping: bool = True,
        check_val_every_n_epoch: int | None = None,
        n_steps_kl_warmup: int | None = None,
        n_epochs_kl_warmup: int | None = 50,
        adversarial_mixing: bool = True,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **kwargs,
    ):
        """Trains the model using amortized variational inference.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset.
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
            Whether to shuffle indices before splitting. If `False`, the val, train, and test set
            are split in the sequential order of the data according to `validation_size` and
            `train_size` percentages.
        batch_size
            Minibatch size to use during training.
        weight_decay
            weight decay regularization term for optimization
        eps
            Optimizer eps
        early_stopping
            Whether to perform early stopping with respect to the validation set.
        check_val_every_n_epoch
            Check val every n train epochs. By default, val is not checked, unless `early_stopping`
            is `True`. If so, val is checked every epoch.
        n_steps_kl_warmup
            Number of training steps (minibatches) to scale weight on KL divergences from 0 to 1.
            Only activated when `n_epochs_kl_warmup` is set to None. If `None`, defaults
            to `floor(0.75 * adata.n_obs)`.
        n_epochs_kl_warmup
            Number of epochs to scale weight on KL divergences from 0 to 1.
            Overrides `n_steps_kl_warmup` when both are not `None`.
        adversarial_mixing
            Whether to use adversarial training to penalize the model for umbalanced mixing of
            modalities.
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
            "adversarial_classifier": adversarial_mixing,
            "weight_decay": weight_decay,
            "eps": eps,
            "n_epochs_kl_warmup": n_epochs_kl_warmup,
            "n_steps_kl_warmup": n_steps_kl_warmup,
            "optimizer": "AdamW",
            "scale_adversarial_loss": 1,
        }
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict

        datasplitter_kwargs = datasplitter_kwargs or {}

        data_splitter = self._data_splitter_cls(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            distributed_sampler=use_distributed_sampler(kwargs.get("strategy", None)),
            batch_size=batch_size,
            **datasplitter_kwargs,
        )
        training_plan = self._training_plan_cls(self.module, **plan_kwargs)
        runner = self._train_runner_cls(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
            early_stopping=early_stopping,
            check_val_every_n_epoch=check_val_every_n_epoch,
            early_stopping_monitor="reconstruction_loss_validation",
            early_stopping_patience=50,
            **kwargs,
        )
        return runner()

    @torch.inference_mode()
    def get_library_size_factors(
        self,
        adata: AnnOrMuData | None = None,
        indices: Sequence[int] = None,
        batch_size: int = 128,
    ) -> dict[str, np.ndarray]:
        """Return library size factors.

        Parameters
        ----------
        adata
            AnnOrMuData object with equivalent structure to initial AnnData. If `None`, defaults
            to the AnnOrMuData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        Library size factor for expression and accessibility
        """
        self._check_adata_modality_weights(adata)
        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        lib_exp = []
        lib_acc = []
        for tensors in scdl:
            outputs = self.module.inference(**self.module._get_inference_input(tensors))
            lib_exp.append(outputs["libsize_expr"].cpu())
            lib_acc.append(outputs["libsize_acc"].cpu())

        return {
            "expression": torch.cat(lib_exp).numpy().squeeze(),
            "accessibility": torch.cat(lib_acc).numpy().squeeze(),
        }

    @torch.inference_mode()
    def get_region_factors(self) -> np.ndarray:
        """Return region-specific factors."""
        if self.n_regions == 0:
            return np.zeros(1)
        else:
            if self.module.region_factors is None:
                raise RuntimeError("region factors were not included in this model")
            return torch.sigmoid(self.module.region_factors).cpu().numpy()

    @torch.inference_mode()
    def get_latent_representation(
        self,
        adata: AnnOrMuData | None = None,
        modality: Literal["joint", "expression", "accessibility"] = "joint",
        indices: Sequence[int] | None = None,
        give_mean: bool = True,
        batch_size: int | None = None,
    ) -> np.ndarray:
        r"""Return the latent representation for each cell.

        Parameters
        ----------
        adata
            AnnOrMuData object with equivalent structure to initial AnnData. If `None`, defaults
            to the AnnOrMuData object used to initialize the model.
        modality
            Return modality specific or joint latent representation.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Give mean of distribution or sample from it.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.

        Returns
        -------
        latent_representation : np.ndarray
            Low-dimensional representation for each cell
        """
        if not self.is_trained_:
            raise RuntimeError("Please train the model first.")
        self._check_adata_modality_weights(adata)
        keys = {"z": "z", "qz_m": "qz_m", "qz_v": "qz_v"}
        if self.fully_paired and modality != "joint":
            raise RuntimeError("A fully paired model only has a joint latent representation.")
        if not self.fully_paired and modality != "joint":
            if modality == "expression":
                keys = {"z": "z_expr", "qz_m": "qzm_expr", "qz_v": "qzv_expr"}
            elif modality == "accessibility":
                keys = {"z": "z_acc", "qz_m": "qzm_acc", "qz_v": "qzv_acc"}
            elif modality == "protein":
                keys = {"z": "z_pro", "qz_m": "qzm_pro", "qz_v": "qzv_pro"}
            else:
                raise RuntimeError(
                    "modality must be 'joint', 'expression', 'accessibility', or 'protein'."
                )

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        latent = []
        for tensors in scdl:
            inference_inputs = self.module._get_inference_input(tensors)
            outputs = self.module.inference(**inference_inputs)
            qz_m = outputs[keys["qz_m"]]
            qz_v = outputs[keys["qz_v"]]
            z = outputs[keys["z"]]

            if give_mean:
                # does each model need to have this latent distribution param?
                if self.module.latent_distribution == "ln":
                    samples = Normal(qz_m, qz_v.sqrt()).sample([1])
                    z = torch.nn.functional.softmax(samples, dim=-1)
                    z = z.mean(dim=0)
                else:
                    z = qz_m

            latent += [z.cpu()]
        return torch.cat(latent).numpy()

    @torch.inference_mode()
    def get_normalized_accessibility(
        self,
        adata: AnnOrMuData | None = None,
        indices: Sequence[int] = None,
        n_samples_overall: int | None = None,
        region_list: Sequence[str] | None = None,
        transform_batch: str | int | None = None,
        use_z_mean: bool = True,
        threshold: float | None = None,
        normalize_cells: bool = False,
        normalize_regions: bool = False,
        batch_size: int = 128,
        return_numpy: bool = False,
    ) -> np.ndarray | csr_matrix | pd.DataFrame:
        """Impute the full accessibility matrix.

        Returns a matrix of accessibility probabilities for each cell and genomic region in the
        input (for return matrix A, A[i,j] is the probability that region j is accessible in cell
        i).

        Parameters
        ----------
        adata
            AnnOrMuData object that has been registered with scvi. If `None`, defaults to the
            AnnOrMuData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples_overall
            Number of samples to return in total
        region_list
            Regions to use. if `None`, all regions are used.
        transform_batch
            Batch to condition on.
            If transform_batch is:

            - None, then real observed batch is used
            - int, then batch transform_batch is used
        use_z_mean
            If True (default), use the distribution mean. Otherwise, sample from the distribution.
        threshold
            If provided, values below the threshold are replaced with 0 and a sparse matrix
            is returned instead. This is recommended for very large matrices. Must be between 0 and
            1.
        normalize_cells
            Whether to reintroduce library size factors to scale the normalized probabilities.
            This makes the estimates closer to the input, but removes the library size correction.
            False by default.
        normalize_regions
            Whether to reintroduce region factors to scale the normalized probabilities. This makes
            the estimates closer to the input, but removes the region-level bias correction. False
            by default.
        batch_size
            Minibatch size for data loading into model
        """
        self._check_adata_modality_weights(adata)
        adata = self._validate_anndata(adata)
        adata_manager = self.get_anndata_manager(adata, required=True)
        if indices is None:
            indices = np.arange(adata.n_obs)
        if n_samples_overall is not None:
            indices = np.random.choice(indices, n_samples_overall)
        post = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        transform_batch = _get_batch_code_from_category(adata_manager, transform_batch)

        if region_list is None:
            region_mask = slice(None)
        else:
            region_mask = [region in region_list for region in adata.var_names[: self.n_regions]]

        if threshold is not None and (threshold < 0 or threshold > 1):
            raise ValueError("the provided threshold must be between 0 and 1")

        imputed = []
        for tensors in post:
            get_generative_input_kwargs = {"transform_batch": transform_batch[0]}
            generative_kwargs = {"use_z_mean": use_z_mean}
            inference_outputs, generative_outputs = self.module.forward(
                tensors=tensors,
                get_generative_input_kwargs=get_generative_input_kwargs,
                generative_kwargs=generative_kwargs,
                compute_loss=False,
            )
            p = generative_outputs["p"].cpu()

            if normalize_cells:
                p *= inference_outputs["libsize_acc"].cpu()
            if normalize_regions:
                p *= torch.sigmoid(self.module.region_factors).cpu()
            if threshold:
                p[p < threshold] = 0
                p = csr_matrix(p.numpy())
            if region_mask is not None:
                p = p[:, region_mask]
            imputed.append(p)

        if threshold:  # imputed is a list of csr_matrix objects
            imputed = vstack(imputed, format="csr")
        else:  # imputed is a list of tensors
            imputed = torch.cat(imputed).numpy()

        if np.all(imputed is None):
            return pd.DataFrame(
                imputed,
                index=adata.obs_names[indices],
                columns=[],
            )
        else:
            if return_numpy:
                return imputed
            elif threshold:
                return pd.DataFrame.sparse.from_spmatrix(
                    imputed,
                    index=adata.obs_names[indices],
                    columns=adata["rna"].var_names[: self.n_regions][region_mask]
                    if isinstance(adata, MuData)
                    else adata.var_names[: self.n_regions][region_mask],
                )
            else:
                return pd.DataFrame(
                    imputed,
                    index=adata.obs_names[indices],
                    columns=adata["rna"].var_names[: self.n_regions][region_mask]
                    if isinstance(adata, MuData)
                    else adata.var_names[: self.n_regions][region_mask],
                )

    @torch.inference_mode()
    def get_normalized_expression(
        self,
        adata: AnnOrMuData | None = None,
        indices: Sequence[int] | None = None,
        n_samples_overall: int | None = None,
        transform_batch: Sequence[Number | str] | None = None,
        gene_list: Sequence[str] | None = None,
        use_z_mean: bool = True,
        n_samples: int = 1,
        batch_size: int | None = None,
        return_mean: bool = True,
        return_numpy: bool = False,
        silent: bool = True,
    ) -> np.ndarray | pd.DataFrame:
        r"""Returns the normalized (decoded) gene expression.

        This is denoted as :math:`\rho_n` in the scVI paper.

        Parameters
        ----------
        adata
            AnnOrMuData object with equivalent structure to initial AnnData. If `None`, defaults
            to the AnnOrMuData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        n_samples_overall
            Number of observations to sample from ``indices`` if ``indices`` is provided.
        transform_batch
            Batch to condition on.
            If transform_batch is:

            - None, then real observed batch is used.
            - int, then batch transform_batch is used.
        gene_list
            Return frequencies of expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        use_z_mean
            If True, use the mean of the latent distribution, otherwise sample from it
        n_samples
            Number of posterior samples to use for estimation.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a numpy array instead of a pandas DataFrame.
        %(de_silent)s

        Returns
        -------
        If `n_samples` > 1 and `return_mean` is False, then the shape is `(samples, cells, genes)`.
        Otherwise, shape is `(cells, genes)`. In this case, return type is
        :class:`~pandas.DataFrame` unless `return_numpy` is True.
        """
        self._check_adata_modality_weights(adata)
        adata = self._validate_anndata(adata)
        adata_manager = self.get_anndata_manager(adata, required=True)
        if indices is None:
            indices = np.arange(adata.n_obs)
        if n_samples_overall is not None:
            indices = np.random.choice(indices, n_samples_overall)
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        transform_batch = _get_batch_code_from_category(adata_manager, transform_batch)

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names[: self.n_genes]
            gene_mask = [gene in gene_list for gene in all_genes]

        exprs = []
        for tensors in scdl:
            per_batch_exprs = []
            for batch in track(transform_batch, disable=silent):
                if batch is not None:
                    batch_indices = tensors[REGISTRY_KEYS.BATCH_KEY]
                    tensors[REGISTRY_KEYS.BATCH_KEY] = torch.ones_like(batch_indices) * batch
                _, generative_outputs = self.module.forward(
                    tensors=tensors,
                    inference_kwargs={"n_samples": n_samples},
                    generative_kwargs={"use_z_mean": use_z_mean},
                    compute_loss=False,
                )
                output = generative_outputs["px_scale"]
                output = output[..., gene_mask]
                output = output.cpu().numpy()
                per_batch_exprs.append(output)
            per_batch_exprs = np.stack(
                per_batch_exprs
            )  # shape is (len(transform_batch) x batch_size x n_var)
            exprs += [per_batch_exprs.mean(0)]

        if n_samples > 1:
            # The -2 axis correspond to cells.
            exprs = np.concatenate(exprs, axis=-2)
        else:
            exprs = np.concatenate(exprs, axis=0)
        if n_samples > 1 and return_mean:
            exprs = exprs.mean(0)

        if return_numpy:
            return exprs
        else:
            return pd.DataFrame(
                exprs,
                columns=adata.var_names[: self.n_genes][gene_mask],
                index=adata.obs_names[indices],
            )

    @de_dsp.dedent
    def differential_accessibility(
        self,
        adata: AnnData | None = None,
        groupby: str | None = None,
        group1: Iterable[str] | None = None,
        group2: str | None = None,
        idx1: Sequence[int] | Sequence[bool] | None = None,
        idx2: Sequence[int] | Sequence[bool] | None = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float = 0.05,
        batch_size: int | None = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Iterable[str] | None = None,
        batchid2: Iterable[str] | None = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        **kwargs,
    ) -> pd.DataFrame:
        r"""A unified method for differential accessibility analysis.

        Implements ``'vanilla'`` DE :cite:p:`Lopez18` and ``'change'`` mode DE :cite:p:`Boyeau19`.

        Parameters
        ----------
        %(de_adata)s
        %(de_groupby)s
        %(de_group1)s
        %(de_group2)s
        %(de_idx1)s
        %(de_idx2)s
        %(de_mode)s
        %(de_delta)s
        %(de_batch_size)s
        %(de_all_stats)s
        %(de_batch_correction)s
        %(de_batchid1)s
        %(de_batchid2)s
        %(de_fdr_target)s
        %(de_silent)s
        **kwargs
            Keyword args for :meth:`scvi.model.base.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential accessibility DataFrame with the following columns:
        prob_da
            the probability of the region being differentially accessible
        is_da_fdr
            whether the region passes a multiple hypothesis correction procedure with the
            target_fdr threshold
        bayes_factor
            Bayes Factor indicating the level of significance of the analysis
        effect_size
            the effect size, computed as (accessibility in population 2) - (accessibility in
            population 1)
        emp_effect
            the empirical effect, based on observed detection rates instead of the estimated
            accessibility scores from the PeakVI model
        est_prob1
            the estimated probability of accessibility in population 1
        est_prob2
            the estimated probability of accessibility in population 2
        emp_prob1
            the empirical (observed) probability of accessibility in population 1
        emp_prob2
            the empirical (observed) probability of accessibility in population 2

        """
        self._check_adata_modality_weights(adata)
        adata = self._validate_anndata(adata)
        col_names = adata.var_names[: self.n_genes]
        model_fn = partial(
            self.get_normalized_accessibility, use_z_mean=False, batch_size=batch_size
        )

        all_stats_fn = partial(
            scatac_raw_counts_properties,
            var_idx=np.arange(adata.shape[1])[: self.n_genes],
        )

        result = _de_core(
            adata_manager=self.get_anndata_manager(adata, required=True),
            model_fn=model_fn,
            representation_fn=None,
            groupby=groupby,
            group1=group1,
            group2=group2,
            idx1=idx1,
            idx2=idx2,
            all_stats=all_stats,
            all_stats_fn=all_stats_fn,
            col_names=col_names,
            mode=mode,
            batchid1=batchid1,
            batchid2=batchid2,
            delta=delta,
            batch_correction=batch_correction,
            fdr=fdr_target,
            silent=silent,
            pseudocounts=1e-6,
            **kwargs,
        )

        # manually change the results DataFrame to fit a PeakVI differential accessibility results
        result = pd.DataFrame(
            {
                "prob_da": result.proba_de,
                "is_da_fdr": result.loc[:, f"is_de_fdr_{fdr_target}"],
                "bayes_factor": result.bayes_factor,
                "effect_size": result.scale2 - result.scale1,
                "emp_effect": result.emp_mean2 - result.emp_mean1,
                "est_prob1": result.scale1,
                "est_prob2": result.scale2,
                "emp_prob1": result.emp_mean1,
                "emp_prob2": result.emp_mean2,
            },
            index=col_names,
        )
        return result

    @de_dsp.dedent
    def differential_expression(
        self,
        adata: AnnData | None = None,
        groupby: str | None = None,
        group1: Iterable[str] | None = None,
        group2: str | None = None,
        idx1: Sequence[int] | Sequence[bool] | None = None,
        idx2: Sequence[int] | Sequence[bool] | None = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float = 0.25,
        batch_size: int | None = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Iterable[str] | None = None,
        batchid2: Iterable[str] | None = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        **kwargs,
    ) -> pd.DataFrame:
        r"""A unified method for differential expression analysis.

        Implements `"vanilla"` DE :cite:p:`Lopez18` and `"change"` mode DE :cite:p:`Boyeau19`.

        Parameters
        ----------
        %(de_adata)s
        %(de_groupby)s
        %(de_group1)s
        %(de_group2)s
        %(de_idx1)s
        %(de_idx2)s
        %(de_mode)s
        %(de_delta)s
        %(de_batch_size)s
        %(de_all_stats)s
        %(de_batch_correction)s
        %(de_batchid1)s
        %(de_batchid2)s
        %(de_fdr_target)s
        %(de_silent)s
        **kwargs
            Keyword args for :meth:`scvi.model.base.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential expression DataFrame.
        """
        self._check_adata_modality_weights(adata)
        adata = self._validate_anndata(adata)

        col_names = adata.var_names[: self.n_genes]
        model_fn = partial(
            self.get_normalized_expression,
            batch_size=batch_size,
        )
        all_stats_fn = partial(
            scrna_raw_counts_properties,
            var_idx=np.arange(adata.shape[1])[: self.n_genes],
        )
        result = _de_core(
            adata_manager=self.get_anndata_manager(adata, required=True),
            model_fn=model_fn,
            representation_fn=None,
            groupby=groupby,
            group1=group1,
            group2=group2,
            idx1=idx1,
            idx2=idx2,
            all_stats=all_stats,
            all_stats_fn=all_stats_fn,
            col_names=col_names,
            mode=mode,
            batchid1=batchid1,
            batchid2=batchid2,
            delta=delta,
            batch_correction=batch_correction,
            fdr=fdr_target,
            silent=silent,
            pseudocounts=1e-6,
            **kwargs,
        )

        return result

    @torch.no_grad()
    def get_protein_foreground_probability(
        self,
        adata: AnnOrMuData | None = None,
        indices: Sequence[int] | None = None,
        transform_batch: Sequence[Number | str] | None = None,
        protein_list: Sequence[str] | None = None,
        n_samples: int = 1,
        batch_size: int | None = None,
        use_z_mean: bool = True,
        return_mean: bool = True,
        return_numpy: bool | None = None,
        silent: bool = True,
    ):
        r"""Returns the foreground probability for proteins.

        This is denoted as :math:`(1 - \pi_{nt})` in the totalVI paper.

        Parameters
        ----------
        adata
            AnnOrMuData object with equivalent structure to initial AnnData. If ``None``, defaults
            to the AnnOrMuData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        transform_batch
            Batch to condition on.
            If transform_batch is:

            * ``None`` - real observed batch is used
            * ``int`` - batch transform_batch is used
            * ``List[int]`` - average over batches in list
        protein_list
            Return protein expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        n_samples
            Number of posterior samples to use for estimation.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame
            includes gene names as columns. If either ``n_samples=1`` or ``return_mean=True``,
            defaults to ``False``. Otherwise, it defaults to `True`.
        %(de_silent)s

        Returns
        -------
        - **foreground_probability** - probability foreground for each protein

        If `n_samples` > 1 and `return_mean` is False, then the shape is `(samples, cells, genes)`.
        Otherwise, shape is `(cells, genes)`. In this case, return type is
        :class:`~pandas.DataFrame` unless `return_numpy` is True.
        """
        adata = self._validate_anndata(adata)
        post = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        if protein_list is None:
            protein_mask = slice(None)
        else:
            all_proteins = self.scvi_setup_dict_["protein_names"]
            protein_mask = [True if p in protein_list else False for p in all_proteins]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                warnings.warn(
                    "`return_numpy` must be `True` if `n_samples > 1` and `return_mean` is "
                    "`False`, returning an `np.ndarray`.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
            return_numpy = True
        if indices is None:
            indices = np.arange(adata.n_obs)

        py_mixings = []
        if not isinstance(transform_batch, IterableClass):
            transform_batch = [transform_batch]

        transform_batch = _get_batch_code_from_category(self.adata_manager, transform_batch)
        for tensors in post:
            y = tensors[REGISTRY_KEYS.PROTEIN_EXP_KEY]
            py_mixing = torch.zeros_like(y[..., protein_mask])
            if n_samples > 1:
                py_mixing = torch.stack(n_samples * [py_mixing])
            for b in track(transform_batch, disable=silent):
                # generative_kwargs = dict(transform_batch=b)
                generative_kwargs = {"use_z_mean": use_z_mean, "transform_batch": b}
                inference_kwargs = {"n_samples": n_samples}
                _, generative_outputs = self.module.forward(
                    tensors=tensors,
                    inference_kwargs=inference_kwargs,
                    generative_kwargs=generative_kwargs,
                    compute_loss=False,
                )
                py_mixing += torch.sigmoid(generative_outputs["py_"]["mixing"])[
                    ..., protein_mask
                ].cpu()
            py_mixing /= len(transform_batch)
            py_mixings += [py_mixing]
        if n_samples > 1:
            # concatenate along batch dimension -> result shape = (samples, cells, features)
            py_mixings = torch.cat(py_mixings, dim=1)
            # (cells, features, samples)
            py_mixings = py_mixings.permute(1, 2, 0)
        else:
            py_mixings = torch.cat(py_mixings, dim=0)

        if return_mean is True and n_samples > 1:
            py_mixings = torch.mean(py_mixings, dim=-1)

        py_mixings = py_mixings.cpu().numpy()

        if return_numpy is True:
            return 1 - py_mixings
        else:
            pro_names = self.protein_state_registry.column_names
            foreground_prob = pd.DataFrame(
                1 - py_mixings,
                columns=pro_names[protein_mask],
                index=adata.obs_names[indices],
            )
            return foreground_prob

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        protein_expression_obsm_key: str | None = None,
        protein_names_uns_key: str | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        %(param_batch_key)s
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        protein_expression_obsm_key
            key in `adata.obsm` for protein expression data.
        protein_names_uns_key
            key in `adata.uns` for protein names. If None, will use the column names of
            `adata.obsm[protein_expression_obsm_key]` if it is a DataFrame, else will assign
            sequential names to proteins.
        """
        warnings.warn(
            "MULTIVI is supposed to work with MuData. the use of anndata is "
            "deprecated and will be removed in scvi-tools 1.4. Please use setup_mudata",
            DeprecationWarning,
            stacklevel=settings.warnings_stacklevel,
        )
        setup_method_args = cls._get_setup_method_args(**locals())
        adata.obs["_indices"] = np.arange(adata.n_obs)
        batch_field = CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key)
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            batch_field,
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, None),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            NumericalJointObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
            NumericalObsField(REGISTRY_KEYS.INDICES_KEY, "_indices"),
        ]
        if protein_expression_obsm_key is not None:
            anndata_fields.append(
                ProteinObsmField(
                    REGISTRY_KEYS.PROTEIN_EXP_KEY,
                    protein_expression_obsm_key,
                    use_batch_mask=True,
                    batch_field=batch_field,
                    colnames_uns_key=protein_names_uns_key,
                    is_count_data=True,
                )
            )

        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    def _check_adata_modality_weights(self, adata):
        """Checks if adata is None and weights are per cell.

        :param adata: anndata object
        :return:
        """
        if (adata is not None) and (self.module.modality_weights == "cell"):
            raise RuntimeError("Held out data not permitted when using per cell weights")

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_mudata(
        cls,
        mdata: MuData,
        rna_layer: str | None = None,
        atac_layer: str | None = None,
        protein_layer: str | None = None,
        batch_key: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        idx_layer: str | None = None,
        modalities: dict[str, str] | None = None,
        **kwargs,
    ):
        """%(summary_mdata)s.

        Parameters
        ----------
        %(param_mdata)s
        rna_layer
            RNA layer key. If `None`, will use `.X` of specified modality key.
        atac_layer
            ATAC layer key. If `None`, will use `.X` of specified modality key.
        protein_layer
            Protein layer key. If `None`, will use `.X` of specified modality key.
        %(param_batch_key)s
        size_factor_key
            Key in `mdata.obsm` for size factors. The first column corresponds to RNA size factors,
            the second to ATAC size factors.
            The second column need to be normalized and between 0 and 1.
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        %(idx_layer)s
        %(param_modalities)s

        Examples
        --------
        >>> mdata = muon.read_10x_h5("filtered_feature_bc_matrix.h5")
        >>> scvi.model.MULTIVI.setup_mudata(
                mdata, modalities={"rna_layer": "rna", "protein_layer": "atac"}
            )
        >>> vae = scvi.model.MULTIVI(mdata)
        """
        setup_method_args = cls._get_setup_method_args(**locals())

        if modalities is None:
            raise ValueError("Modalities cannot be None.")
        modalities = cls._create_modalities_attr_dict(modalities, setup_method_args)
        mdata.obs["_indices"] = np.arange(mdata.n_obs)

        batch_field = fields.MuDataCategoricalObsField(
            REGISTRY_KEYS.BATCH_KEY,
            batch_key,
            mod_key=modalities.batch_key,
        )
        mudata_fields = [
            batch_field,
            fields.MuDataCategoricalObsField(
                REGISTRY_KEYS.LABELS_KEY,
                None,
                mod_key=None,
            ),
            fields.MuDataNumericalJointObsField(
                REGISTRY_KEYS.SIZE_FACTOR_KEY,
                size_factor_key,
                mod_key=None,
                required=False,
            ),
            fields.MuDataCategoricalJointObsField(
                REGISTRY_KEYS.CAT_COVS_KEY,
                categorical_covariate_keys,
                mod_key=modalities.categorical_covariate_keys,
            ),
            fields.MuDataNumericalJointObsField(
                REGISTRY_KEYS.CONT_COVS_KEY,
                continuous_covariate_keys,
                mod_key=modalities.continuous_covariate_keys,
            ),
            fields.MuDataNumericalObsField(
                REGISTRY_KEYS.INDICES_KEY,
                "_indices",
                mod_key=modalities.idx_layer,
                required=False,
            ),
        ]
        if modalities.rna_layer is not None:
            mudata_fields.append(
                fields.MuDataLayerField(
                    REGISTRY_KEYS.X_KEY,
                    rna_layer,
                    mod_key=modalities.rna_layer,
                    is_count_data=True,
                    mod_required=True,
                )
            )
        if modalities.atac_layer is not None:
            mudata_fields.append(
                fields.MuDataLayerField(
                    REGISTRY_KEYS.ATAC_X_KEY,
                    atac_layer,
                    mod_key=modalities.atac_layer,
                    is_count_data=True,
                    mod_required=True,
                )
            )
        if modalities.protein_layer is not None:
            mudata_fields.append(
                fields.MuDataProteinLayerField(
                    REGISTRY_KEYS.PROTEIN_EXP_KEY,
                    protein_layer,
                    mod_key=modalities.protein_layer,
                    use_batch_mask=True,
                    batch_field=batch_field,
                    is_count_data=True,
                    mod_required=True,
                )
            )
        adata_manager = AnnDataManager(fields=mudata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(mdata, **kwargs)
        cls.register_manager(adata_manager)
