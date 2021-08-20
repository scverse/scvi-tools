import logging
from functools import partial
from typing import Iterable, Optional, Sequence, Union

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from scipy.sparse import csr_matrix, vstack
from torch.distributions import Normal

from scvi import _CONSTANTS
from scvi._compat import Literal
from scvi._docs import doc_differential_expression
from scvi._utils import _doc_params
from scvi.dataloaders import DataSplitter
from scvi.model._utils import (
    _get_batch_code_from_category,
    _get_var_names_from_setup_anndata,
    scatac_raw_counts_properties,
    scrna_raw_counts_properties,
)
from scvi.model.base import UnsupervisedTrainingMixin
from scvi.module import MULTIVAE
from scvi.train import AdversarialTrainingPlan, TrainRunner
from scvi.train._callbacks import SaveBestState

from .base import BaseModelClass, VAEMixin
from .base._utils import _de_core

Number = Union[int, float]
logger = logging.getLogger(__name__)


class MULTIVI(VAEMixin, UnsupervisedTrainingMixin, BaseModelClass):
    """
    Integration of multi-modal and single-modality data [AshuachGabitto21]_.

    MultiVI is used to integrate multiomic datasets with single-modality (expression
    or accessibility) datasets.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    n_genes
        The number of gene expression features (genes).
    n_regions
        The number of accessibility features (genomic regions).
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
        Keyword args for :class:`~scvi.module.PEAKVAE`

    Examples
    --------
    >>> adata_rna = anndata.read_h5ad(path_to_rna_anndata)
    >>> adata_atac = scvi.data.read_10x_atac(path_to_atac_anndata)
    >>> adata_multi = scvi.data.read_10x_multiome(path_to_multiomic_anndata)
    >>> adata_mvi = scvi.data.organize_multiome_anndatas(adata_multi, adata_rna, adata_atac)
    >>> scvi.data.setup_anndata(adata_mvi, batch_key="modality")
    >>> vae = scvi.model.MULTIVI(adata_mvi)
    >>> vae.train()

    Notes
    ------
    * The model assumes that the features are organized so that all expression features are
       consecutive, followed by all accessibility features. For example, if the data has 100 genes
       and 250 genomic regions, the model assumes that the first 100 features are genes, and the
       next 250 are the regions.

    * The main batch annotation, specified in the `scvi.data.setup_anndata`, should correspond to
       the modality each cell originated from. This allows the model to focus mixing efforts, using
       an adversarial component, on mixing the modalities. Other covariates can be specified using
       the `categorical_covariate_keys` argument.
    """

    def __init__(
        self,
        adata: AnnData,
        n_genes: int,
        n_regions: int,
        n_hidden: Optional[int] = None,
        n_latent: Optional[int] = None,
        n_layers_encoder: int = 2,
        n_layers_decoder: int = 2,
        dropout_rate: float = 0.1,
        region_factors: bool = True,
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb",
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        latent_distribution: Literal["normal", "ln"] = "normal",
        deeply_inject_covariates: bool = False,
        encode_covariates: bool = False,
        fully_paired: bool = False,
        **model_kwargs,
    ):
        super().__init__(adata)

        n_cats_per_cov = (
            self.scvi_setup_dict_["extra_categoricals"]["n_cats_per_key"]
            if "extra_categoricals" in self.scvi_setup_dict_
            else []
        )

        self.module = MULTIVAE(
            n_input_genes=n_genes,
            n_input_regions=n_regions,
            n_batch=self.summary_stats["n_batch"],
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers_encoder=n_layers_encoder,
            n_layers_decoder=n_layers_decoder,
            n_continuous_cov=self.summary_stats["n_continuous_covs"],
            n_cats_per_cov=n_cats_per_cov,
            dropout_rate=dropout_rate,
            region_factors=region_factors,
            gene_likelihood=gene_likelihood,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            latent_distribution=latent_distribution,
            deeply_inject_covariates=deeply_inject_covariates,
            encode_covariates=encode_covariates,
            **model_kwargs,
        )
        self._model_summary_string = (
            "MultiVI Model with INPUTS: n_genes:{}, n_regions:{}\n"
            "n_hidden: {}, n_latent: {}, n_layers_encoder: {}, "
            "n_layers_decoder: {} , dropout_rate: {}, latent_distribution: {}, deep injection: {}, "
            "gene_likelihood: {}"
        ).format(
            n_genes,
            n_regions,
            self.module.n_hidden,
            self.module.n_latent,
            n_layers_encoder,
            n_layers_decoder,
            dropout_rate,
            latent_distribution,
            deeply_inject_covariates,
            gene_likelihood,
        )
        self.fully_paired = fully_paired
        self.n_latent = n_latent
        self.init_params_ = self._get_init_params(locals())
        self.n_genes = n_genes

    def train(
        self,
        max_epochs: int = 500,
        lr: float = 1e-4,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        weight_decay: float = 1e-3,
        eps: float = 1e-08,
        early_stopping: bool = True,
        save_best: bool = True,
        check_val_every_n_epoch: Optional[int] = None,
        n_steps_kl_warmup: Optional[int] = None,
        n_epochs_kl_warmup: Optional[int] = 50,
        adversarial_mixing: bool = True,
        plan_kwargs: Optional[dict] = None,
        **kwargs,
    ):
        """
        Trains the model using amortized variational inference.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset.
        lr
            Learning rate for optimization.
        use_gpu
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str), or use CPU (if False).
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        weight_decay
            weight decay regularization term for optimization
        eps
            Optimizer eps
        early_stopping
            Whether to perform early stopping with respect to the validation set.
        save_best
            Save the best model state with respect to the validation loss, or use the final
            state in the training procedure
        check_val_every_n_epoch
            Check val every n train epochs. By default, val is not checked, unless `early_stopping` is `True`.
            If so, val is checked every epoch.
        n_steps_kl_warmup
            Number of training steps (minibatches) to scale weight on KL divergences from 0 to 1.
            Only activated when `n_epochs_kl_warmup` is set to None. If `None`, defaults
            to `floor(0.75 * adata.n_obs)`.
        n_epochs_kl_warmup
            Number of epochs to scale weight on KL divergences from 0 to 1.
            Overrides `n_steps_kl_warmup` when both are not `None`.
        adversarial_mixing
            Whether to use adversarial training to penalize the model for umbalanced mixing of modalities.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        update_dict = dict(
            lr=lr,
            adversarial_classifier=adversarial_mixing,
            weight_decay=weight_decay,
            eps=eps,
            n_epochs_kl_warmup=n_epochs_kl_warmup,
            n_steps_kl_warmup=n_steps_kl_warmup,
            check_val_every_n_epoch=check_val_every_n_epoch,
            early_stopping=early_stopping,
            early_stopping_monitor="reconstruction_loss_validation",
            early_stopping_patience=50,
            optimizer="AdamW",
            scale_adversarial_loss=1,
        )
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict

        if save_best:
            if "callbacks" not in kwargs.keys():
                kwargs["callbacks"] = []
            kwargs["callbacks"].append(
                SaveBestState(monitor="reconstruction_loss_validation")
            )

        data_splitter = DataSplitter(
            self.adata,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            use_gpu=use_gpu,
        )
        training_plan = AdversarialTrainingPlan(self.module, **plan_kwargs)
        runner = TrainRunner(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            early_stopping=early_stopping,
            **kwargs,
        )
        return runner()

    @torch.no_grad()
    def get_library_size_factors(
        self,
        adata: Optional[AnnData] = None,
        indices: Sequence[int] = None,
        batch_size: int = 128,
    ):
        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )

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

    @torch.no_grad()
    def get_region_factors(self) -> np.ndarray:
        if self.module.region_factors is None:
            raise RuntimeError("region factors were not included in this model")
        return torch.sigmoid(self.module.region_factors).cpu().numpy()

    @torch.no_grad()
    def get_latent_representation(
        self,
        adata: Optional[AnnData] = None,
        modality: Literal["joint", "expression", "accessibility"] = "joint",
        indices: Optional[Sequence[int]] = None,
        give_mean: bool = True,
        batch_size: Optional[int] = None,
    ) -> np.ndarray:
        if not self.is_trained_:
            raise RuntimeError("Please train the model first.")

        keys = {"z": "z", "qz_m": "qz_m", "qz_v": "qz_v"}
        if self.fully_paired and modality != "joint":
            raise RuntimeError(
                "A fully paired model only has a joint latent representation."
            )
        if not self.fully_paired and modality != "joint":
            if modality == "expression":
                keys = {"z": "z_expr", "qz_m": "qzm_expr", "qz_v": "qzv_expr"}
            elif modality == "accessibility":
                keys = {"z": "z_acc", "qz_m": "qzm_acc", "qz_v": "qzv_acc"}
            else:
                raise RuntimeError(
                    "modality must be 'joint', 'expression', or 'accessibility'."
                )

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )
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

    @torch.no_grad()
    def get_accessibility_estimates(
        self,
        adata: Optional[AnnData] = None,
        indices: Sequence[int] = None,
        n_samples_overall: Optional[int] = None,
        region_indices: Sequence[int] = None,
        transform_batch: Optional[Union[str, int]] = None,
        use_z_mean: bool = True,
        threshold: Optional[float] = None,
        normalize_cells: bool = False,
        normalize_regions: bool = False,
        batch_size: int = 128,
    ) -> Union[np.ndarray, csr_matrix]:
        adata = self._validate_anndata(adata)
        if indices is None:
            indices = np.arange(adata.n_obs)
        if n_samples_overall is not None:
            indices = np.random.choice(indices, n_samples_overall)
        post = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )
        transform_batch = _get_batch_code_from_category(adata, transform_batch)

        if threshold is not None and (threshold < 0 or threshold > 1):
            raise ValueError("the provided threshold must be between 0 and 1")

        imputed = []
        for tensors in post:
            get_generative_input_kwargs = dict(transform_batch=transform_batch[0])
            generative_kwargs = dict(use_z_mean=use_z_mean)
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
            if region_indices is not None:
                p = p[:, region_indices]
            imputed.append(p)

        if threshold:  # imputed is a list of csr_matrix objects
            imputed = vstack(imputed, format="csr")
        else:  # imputed is a list of tensors
            imputed = torch.cat(imputed).numpy()

        return imputed

    @torch.no_grad()
    def get_normalized_expression(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        n_samples_overall: Optional[int] = None,
        transform_batch: Optional[Sequence[Union[Number, str]]] = None,
        gene_list: Optional[Sequence[str]] = None,
        use_z_mean: bool = True,
        n_samples: int = 1,
        batch_size: Optional[int] = None,
        return_mean: bool = True,
    ) -> Union[np.ndarray, pd.DataFrame]:
        r"""
        Returns the normalized (decoded) gene expression.

        This is denoted as :math:`\rho_n` in the scVI paper.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        transform_batch
            Batch to condition on.
            If transform_batch is:

            - None, then real observed batch is used.
            - int, then batch transform_batch is used.
        gene_list
            Return frequencies of expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        library_size
            Scale the expression frequencies to a common library size.
            This allows gene expression levels to be interpreted on a common scale of relevant
            magnitude. If set to `"latent"`, use the latent libary size.
        use_z_mean
            If True, use the mean of the latent distribution, otherwise sample from it
        n_samples
            Number of posterior samples to use for estimation.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.

        Returns
        -------
        If `n_samples` > 1 and `return_mean` is False, then the shape is `(samples, cells, genes)`.
        Otherwise, shape is `(cells, genes)`. In this case, return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        """
        adata = self._validate_anndata(adata)
        if indices is None:
            indices = np.arange(adata.n_obs)
        if n_samples_overall is not None:
            indices = np.random.choice(indices, n_samples_overall)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )

        transform_batch = _get_batch_code_from_category(adata, transform_batch)

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = _get_var_names_from_setup_anndata(adata)
            gene_mask = [gene in gene_list for gene in all_genes]

        exprs = []
        for tensors in scdl:
            per_batch_exprs = []
            for batch in transform_batch:
                if batch is not None:
                    batch_indices = tensors[_CONSTANTS.BATCH_KEY]
                    tensors[_CONSTANTS.BATCH_KEY] = (
                        torch.ones_like(batch_indices) * batch
                    )
                _, generative_outputs = self.module.forward(
                    tensors=tensors,
                    inference_kwargs=dict(n_samples=n_samples),
                    generative_kwargs=dict(use_z_mean=use_z_mean),
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
        return exprs

    @_doc_params(doc_differential_expression=doc_differential_expression)
    def differential_accessibility(
        self,
        adata: Optional[AnnData] = None,
        groupby: Optional[str] = None,
        group1: Optional[Iterable[str]] = None,
        group2: Optional[str] = None,
        idx1: Optional[Union[Sequence[int], Sequence[bool]]] = None,
        idx2: Optional[Union[Sequence[int], Sequence[bool]]] = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float = 0.05,
        batch_size: Optional[int] = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Optional[Iterable[str]] = None,
        batchid2: Optional[Iterable[str]] = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        two_sided: bool = True,
        **kwargs,
    ) -> pd.DataFrame:
        r"""
        A unified method for differential accessibility analysis.

        Implements `"vanilla"` DE [Lopez18]_ and `"change"` mode DE [Boyeau19]_.

        Parameters
        ----------
        {doc_differential_expression}
        two_sided
            Whether to perform a two-sided test, or a one-sided test.
        **kwargs
            Keyword args for :func:`scvi.utils.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential accessibility DataFrame with the following columns:
        prob_da
            the probability of the region being differentially accessible
        is_da_fdr
            whether the region passes a multiple hypothesis correction procedure with the target_fdr
            threshold
        bayes_factor
            Bayes Factor indicating the level of significance of the analysis
        effect_size
            the effect size, computed as (accessibility in population 2) - (accessibility in population 1)
        emp_effect
            the empirical effect, based on observed detection rates instead of the estimated accessibility
            scores from the PeakVI model
        est_prob1
            the estimated probability of accessibility in population 1
        est_prob2
            the estimated probability of accessibility in population 2
        emp_prob1
            the empirical (observed) probability of accessibility in population 1
        emp_prob2
            the empirical (observed) probability of accessibility in population 2

        """
        adata = self._validate_anndata(adata)
        col_names = _get_var_names_from_setup_anndata(adata)[self.n_genes :]
        model_fn = partial(
            self.get_accessibility_estimates, use_z_mean=False, batch_size=batch_size
        )

        # TODO check if change_fn in kwargs and raise error if so
        def change_fn(a, b):
            return a - b

        if two_sided:

            def m1_domain_fn(samples):
                return np.abs(samples) >= delta

        else:

            def m1_domain_fn(samples):
                return samples >= delta

        all_stats_fn = partial(
            scatac_raw_counts_properties,
            var_idx=np.arange(adata.shape[1])[self.n_genes :],
        )

        result = _de_core(
            adata=adata,
            model_fn=model_fn,
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
            change_fn=change_fn,
            m1_domain_fn=m1_domain_fn,
            silent=silent,
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
        )
        return result

    @_doc_params(doc_differential_expression=doc_differential_expression)
    def differential_expression(
        self,
        adata: Optional[AnnData] = None,
        groupby: Optional[str] = None,
        group1: Optional[Iterable[str]] = None,
        group2: Optional[str] = None,
        idx1: Optional[Union[Sequence[int], Sequence[bool]]] = None,
        idx2: Optional[Union[Sequence[int], Sequence[bool]]] = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float = 0.25,
        batch_size: Optional[int] = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Optional[Iterable[str]] = None,
        batchid2: Optional[Iterable[str]] = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        **kwargs,
    ) -> pd.DataFrame:
        r"""
        A unified method for differential expression analysis.

        Implements `"vanilla"` DE [Lopez18]_ and `"change"` mode DE [Boyeau19]_.

        Parameters
        ----------
        {doc_differential_expression}
        **kwargs
            Keyword args for :func:`scvi.utils.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential expression DataFrame.
        """
        adata = self._validate_anndata(adata)

        col_names = _get_var_names_from_setup_anndata(adata)[: self.n_genes]
        model_fn = partial(
            self.get_normalized_expression,
            batch_size=batch_size,
        )
        all_stats_fn = partial(
            scrna_raw_counts_properties,
            var_idx=np.arange(adata.shape[1])[: self.n_genes],
        )
        result = _de_core(
            adata,
            model_fn=model_fn,
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
            **kwargs,
        )

        return result
