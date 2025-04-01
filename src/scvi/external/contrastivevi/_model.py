"""Model class for contrastive-VI for single cell expression data."""

from __future__ import annotations

import logging
import warnings
from functools import partial
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import torch

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
)
from scvi.dataloaders import AnnDataLoader
from scvi.model._utils import (
    _get_batch_code_from_category,
    _init_library_size,
    get_max_epochs_heuristic,
    scrna_raw_counts_properties,
    use_distributed_sampler,
)
from scvi.model.base import BaseModelClass
from scvi.model.base._de_core import _de_core
from scvi.train import TrainingPlan, TrainRunner
from scvi.utils import setup_anndata_dsp, track
from scvi.utils._docstrings import devices_dsp

from ._contrastive_data_splitting import ContrastiveDataSplitter
from ._module import ContrastiveVAE

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from anndata import AnnData

logger = logging.getLogger(__name__)
Number = int | float


class ContrastiveVI(BaseModelClass):
    """contrastive variational inference :cite:p:`Weinberger23`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via
        :meth:`~scvi.model.ContrastiveVI.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_background_latent
        Dimensionality of the background shared latent space.
    n_salient_latent
        Dimensionality of the salient latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    use_observed_lib_size
        Use observed library size for RNA as scaling factor in mean of conditional distribution.
    wasserstein_penalty
        Weight of the Wasserstein distance loss that further discourages background
        shared variations from leaking into the salient latent space.

    Notes
    -----
    See further usage examples in the following tutorial:

    1. :doc:`/tutorials/notebooks/scrna/contrastiveVI_tutorial`
    """

    _module_cls = ContrastiveVAE
    _data_splitter_cls = ContrastiveDataSplitter
    _training_plan_cls = TrainingPlan
    _train_runner_cls = TrainRunner

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_background_latent: int = 10,
        n_salient_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        use_observed_lib_size: bool = True,
        wasserstein_penalty: float = 0,
    ) -> None:
        super().__init__(adata)

        n_cats_per_cov = (
            self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )
        n_batch = self.summary_stats.n_batch

        library_log_means, library_log_vars = None, None
        if not use_observed_lib_size:
            library_log_means, library_log_vars = _init_library_size(self.adata_manager, n_batch)

        self.module = self._module_cls(
            n_input=self.summary_stats.n_vars,
            n_batch=n_batch,
            n_hidden=n_hidden,
            n_background_latent=n_background_latent,
            n_salient_latent=n_salient_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            use_observed_lib_size=use_observed_lib_size,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            wasserstein_penalty=wasserstein_penalty,
        )
        self._model_summary_string = (
            f"ContrastiveVI Model with the following params: \nn_hidden: {n_hidden}, "
            f"n_background_latent: {n_background_latent}, n_salient_latent: {n_salient_latent}, "
            f"n_layers: {n_layers}, dropout_rate: {dropout_rate}, "
            f"use_observed_lib_size: {use_observed_lib_size}, "
            f"wasserstein_penalty: {wasserstein_penalty}"
        )
        self.init_params_ = self._get_init_params(locals())

    @devices_dsp.dedent
    def train(
        self,
        background_indices: list[int],
        target_indices: list[int],
        max_epochs: int | None = None,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float | None = None,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        load_sparse_tensor: bool = False,
        batch_size: int = 128,
        early_stopping: bool = False,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **trainer_kwargs,
    ):
        """Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset. If `None`, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`
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
        load_sparse_tensor
            ``EXPERIMENTAL`` If ``True``, loads data with sparse CSR or CSC layout as a
            :class:`~torch.Tensor` with the same layout. Can lead to speedups in data transfers to
            GPUs, depending on the sparsity of the data.
        batch_size
            Minibatch size to use during training.
        early_stopping
            Perform early stopping. Additional arguments can be passed in `**kwargs`.
            See :class:`~scvi.train.Trainer` for further options.
        datasplitter_kwargs
            Additional keyword arguments passed into
            :class:`~scvi.dataloaders.ContrastiveDataSplitter`.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        if max_epochs is None:
            max_epochs = get_max_epochs_heuristic(self.adata.n_obs)

        plan_kwargs = plan_kwargs or {}
        datasplitter_kwargs = datasplitter_kwargs or {}

        data_splitter = self._data_splitter_cls(
            self.adata_manager,
            background_indices=background_indices,
            target_indices=target_indices,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            shuffle_set_split=shuffle_set_split,
            distributed_sampler=use_distributed_sampler(trainer_kwargs.get("strategy", None)),
            load_sparse_tensor=load_sparse_tensor,
            **datasplitter_kwargs,
        )
        training_plan = self._training_plan_cls(self.module, **plan_kwargs)

        es = "early_stopping"
        trainer_kwargs[es] = (
            early_stopping if es not in trainer_kwargs.keys() else trainer_kwargs[es]
        )
        runner = self._train_runner_cls(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
            **trainer_kwargs,
        )
        return runner()

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        labels_key: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    @torch.inference_mode()
    def get_latent_representation(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        give_mean: bool = True,
        batch_size: int | None = None,
        representation_kind: str = "salient",
    ) -> np.ndarray:
        """Returns the background or salient latent representation for each cell.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`,
            defaults to the AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Give mean of distribution or sample from it.
        batch_size
            Mini-batch size for data loading into model. Defaults to
            `scvi.settings.batch_size`.
        representation_kind
            Either "background" or "salient" for the corresponding representation kind.

        Returns
        -------
            A numpy array with shape `(n_cells, n_latent)`.
        """
        available_representation_kinds = ["background", "salient"]
        if representation_kind not in available_representation_kinds:
            raise ValueError(
                f"representation_kind = {representation_kind} is not one of"
                f" {available_representation_kinds}"
            )

        adata = self._validate_anndata(adata)
        data_loader = self._make_data_loader(
            adata=adata,
            indices=indices,
            batch_size=batch_size,
            shuffle=False,
            data_loader_class=AnnDataLoader,
        )
        latent = []
        for tensors in data_loader:
            x = tensors[REGISTRY_KEYS.X_KEY]
            batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
            outputs = self.module._generic_inference(x=x, batch_index=batch_index, n_samples=1)

            if representation_kind == "background":
                latent_m = outputs["qz_m"]
                latent_sample = outputs["z"]
            else:
                latent_m = outputs["qs_m"]
                latent_sample = outputs["s"]

            if give_mean:
                latent_sample = latent_m

            latent += [latent_sample.detach().cpu()]
        return torch.cat(latent).numpy()

    @torch.inference_mode()
    def get_normalized_expression(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        transform_batch: Sequence[Number | str] | None = None,
        gene_list: Sequence[str] | None = None,
        library_size: float | str = 1.0,
        n_samples: int = 1,
        n_samples_overall: int | None = None,
        batch_size: int | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
        silent: bool = True,
    ) -> dict[str, np.ndarray | pd.DataFrame]:
        """Returns the normalized (decoded) gene expression.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`,
            defaults to the AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        transform_batch
            Batch to condition on. If transform_batch is:
            - None, then real observed batch is used.
            - int, then batch transform_batch is used.
        gene_list
            Return frequencies of expression for a subset of genes. This can save
            memory when working with large datasets and few genes are of interest.
        library_size
            Scale the expression frequencies to a common library size. This
            allows gene expression levels to be interpreted on a common scale of
            relevant magnitude. If set to `"latent"`, use the latent library size.
        n_samples
            Number of posterior samples to use for estimation.
        n_samples_overall
            The number of random samples in `adata` to use.
        batch_size
            Mini-batch size for data loading into model. Defaults to
            `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a `numpy.ndarray` instead of a `pandas.DataFrame`.
            DataFrame includes gene names as columns. If either `n_samples=1` or
            `return_mean=True`, defaults to `False`. Otherwise, it defaults to `True`.
        %(de_silent)s

        Returns
        -------
            A dictionary with keys "background" and "salient", with value as follows.
            If `n_samples` > 1 and `return_mean` is `False`, then the shape is
            `(samples, cells, genes)`. Otherwise, shape is `(cells, genes)`. In this
            case, return type is `pandas.DataFrame` unless `return_numpy` is `True`.
        """
        adata = self._validate_anndata(adata)
        if indices is None:
            indices = np.arange(adata.n_obs)
        if n_samples_overall is not None:
            indices = np.random.choice(indices, n_samples_overall)
        data_loader = self._make_data_loader(
            adata=adata,
            indices=indices,
            batch_size=batch_size,
            shuffle=False,
            data_loader_class=AnnDataLoader,
        )

        transform_batch = _get_batch_code_from_category(
            self.get_anndata_manager(adata, required=True), transform_batch
        )

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                warnings.warn(
                    "return_numpy must be True if n_samples > 1 and"
                    " return_mean is False, returning np.ndarray",
                    stacklevel=settings.warnings_stacklevel,
                )
            return_numpy = True
        if library_size == "latent":
            generative_output_key = "px_rate"
            scaling = 1
        else:
            generative_output_key = "px_scale"
            scaling = library_size

        background_exprs = []
        salient_exprs = []
        for tensors in data_loader:
            x = tensors[REGISTRY_KEYS.X_KEY]
            batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
            background_per_batch_exprs = []
            salient_per_batch_exprs = []
            for batch in track(transform_batch, disable=silent):
                if batch is not None:
                    batch_index = torch.ones_like(batch_index) * batch
                inference_outputs = self.module._generic_inference(
                    x=x, batch_index=batch_index, n_samples=n_samples
                )
                z = inference_outputs["z"]
                s = inference_outputs["s"]
                library = inference_outputs["library"]
                background_generative_outputs = self.module._generic_generative(
                    z=z, s=torch.zeros_like(s), library=library, batch_index=batch_index
                )
                salient_generative_outputs = self.module._generic_generative(
                    z=z, s=s, library=library, batch_index=batch_index
                )
                background_outputs = self._preprocess_normalized_expression(
                    background_generative_outputs,
                    generative_output_key,
                    gene_mask,
                    scaling,
                )
                background_per_batch_exprs.append(background_outputs)
                salient_outputs = self._preprocess_normalized_expression(
                    salient_generative_outputs,
                    generative_output_key,
                    gene_mask,
                    scaling,
                )
                salient_per_batch_exprs.append(salient_outputs)

            background_per_batch_exprs = np.stack(
                background_per_batch_exprs
            )  # Shape is (len(transform_batch) x batch_size x n_var).
            salient_per_batch_exprs = np.stack(salient_per_batch_exprs)
            background_exprs += [background_per_batch_exprs.mean(0)]
            salient_exprs += [salient_per_batch_exprs.mean(0)]

        if n_samples > 1:
            # The -2 axis correspond to cells.
            background_exprs = np.concatenate(background_exprs, axis=-2)
            salient_exprs = np.concatenate(salient_exprs, axis=-2)
        else:
            background_exprs = np.concatenate(background_exprs, axis=0)
            salient_exprs = np.concatenate(salient_exprs, axis=0)
        if n_samples > 1 and return_mean:
            background_exprs = background_exprs.mean(0)
            salient_exprs = salient_exprs.mean(0)

        if return_numpy is None or return_numpy is False:
            genes = adata.var_names[gene_mask]
            samples = adata.obs_names[indices]
            background_exprs = pd.DataFrame(background_exprs, columns=genes, index=samples)
            salient_exprs = pd.DataFrame(salient_exprs, columns=genes, index=samples)
        return {"background": background_exprs, "salient": salient_exprs}

    @torch.inference_mode()
    def get_salient_normalized_expression(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        transform_batch: Sequence[Number | str] | None = None,
        gene_list: Sequence[str] | None = None,
        library_size: float | str = 1.0,
        n_samples: int = 1,
        n_samples_overall: int | None = None,
        batch_size: int | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
    ) -> np.ndarray | pd.DataFrame:
        """Returns the normalized (decoded) gene expression.

        Gene expressions are decoded from both the background and salient latent space.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`,
            defaults to the AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        transform_batch
            Batch to condition on. If transform_batch is:
            - None, then real observed batch is used.
            - int, then batch transform_batch is used.
        gene_list
            Return frequencies of expression for a subset of genes. This can
            save memory when working with large datasets and few genes are of interest.
        library_size
            Scale the expression frequencies to a common library size. This
            allows gene expression levels to be interpreted on a common scale of
            relevant magnitude. If set to `"latent"`, use the latent library size.
        n_samples
            Number of posterior samples to use for estimation.
        n_samples_overall
            The number of random samples in `adata` to use.
        batch_size
            Mini-batch size for data loading into model. Defaults to
            `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a `numpy.ndarray` instead of a `pandas.DataFrame`.
            DataFrame includes gene names as columns. If either `n_samples=1` or
            `return_mean=True`, defaults to `False`. Otherwise, it defaults to `True`.

        Returns
        -------
            If `n_samples` > 1 and `return_mean` is `False`, then the shape is
            `(samples, cells, genes)`. Otherwise, shape is `(cells, genes)`. In this
            case, return type is `pandas.DataFrame` unless `return_numpy` is `True`.
        """
        exprs = self.get_normalized_expression(
            adata=adata,
            indices=indices,
            transform_batch=transform_batch,
            gene_list=gene_list,
            library_size=library_size,
            n_samples=n_samples,
            n_samples_overall=n_samples_overall,
            batch_size=batch_size,
            return_mean=return_mean,
            return_numpy=return_numpy,
        )
        return exprs["salient"]

    @torch.inference_mode()
    def get_specific_normalized_expression(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        transform_batch: Sequence[Number | str] | None = None,
        gene_list: Sequence[str] | None = None,
        library_size: float | str = 1,
        n_samples: int = 1,
        n_samples_overall: int | None = None,
        batch_size: int | None = None,
        return_mean: bool = True,
        return_numpy: bool | None = None,
        expression_type: str | None = None,
        indices_to_return_salient: Sequence[int] | None = None,
    ):
        """Returns the normalized (decoded) gene expression.

        Gene expressions are decoded from either the background or salient latent space.
        One of `expression_type` or `indices_to_return_salient` should have an input
        argument.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`,
            defaults to the AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        transform_batch
            Batch to condition on. If transform_batch is:
            - None, then real observed batch is used.
            - int, then batch transform_batch is used.
        gene_list
            Return frequencies of expression for a subset of genes. This can
            save memory when working with large datasets and few genes are of interest.
        library_size
            Scale the expression frequencies to a common library size. This
            allows gene expression levels to be interpreted on a common scale of
            relevant magnitude. If set to `"latent"`, use the latent library size.
        n_samples
            Number of posterior samples to use for estimation.
        n_samples_overall
            The number of random samples in `adata` to use.
        batch_size
            Mini-batch size for data loading into model. Defaults to
            `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a `numpy.ndarray` instead of a `pandas.DataFrame`.
            DataFrame includes gene names as columns. If either `n_samples=1` or
            `return_mean=True`, defaults to `False`. Otherwise, it defaults to `True`.
        expression_type
            One of {"salient", "background"} to specify the type of
            normalized expression to return.
        indices_to_return_salient
            If `indices` is a subset of `indices_to_return_salient`, normalized
            expressions derived from background and salient latent embeddings are
            returned. If `indices` is not `None` and is not a subset of
            `indices_to_return_salient`, normalized expressions derived only from
            background latent embeddings are returned.

        Returns
        -------
            If `n_samples` > 1 and `return_mean` is `False`, then the shape is
            `(samples, cells, genes)`. Otherwise, shape is `(cells, genes)`. In this
            case, return type is `pandas.DataFrame` unless `return_numpy` is `True`.
        """
        is_expression_type_none = expression_type is None
        is_indices_to_return_salient_none = indices_to_return_salient is None
        if is_expression_type_none and is_indices_to_return_salient_none:
            raise ValueError(
                "Both expression_type and indices_to_return_salient are None! "
                "Exactly one of them needs to be supplied with an input argument."
            )
        elif (not is_expression_type_none) and (not is_indices_to_return_salient_none):
            raise ValueError(
                "Both expression_type and indices_to_return_salient have an input "
                "argument! Exactly one of them needs to be supplied with an input "
                "argument."
            )
        else:
            exprs = self.get_normalized_expression(
                adata=adata,
                indices=indices,
                transform_batch=transform_batch,
                gene_list=gene_list,
                library_size=library_size,
                n_samples=n_samples,
                n_samples_overall=n_samples_overall,
                batch_size=batch_size,
                return_mean=return_mean,
                return_numpy=return_numpy,
            )
            if not is_expression_type_none:
                return exprs[expression_type]
            else:
                if indices is None:
                    indices = np.arange(adata.n_obs)
                if set(indices).issubset(set(indices_to_return_salient)):
                    return exprs["salient"]
                else:
                    return exprs["background"]

    def differential_expression(
        self,
        adata: AnnData | None = None,
        groupby: str | None = None,
        group1: Iterable[str] | None = None,
        group2: str | None = None,
        idx1: Sequence[int] | (Sequence[bool] | str) | None = None,
        idx2: Sequence[int] | (Sequence[bool] | str) | None = None,
        mode: str = "change",
        delta: float = 0.25,
        batch_size: int | None = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Iterable[str] | None = None,
        batchid2: Iterable[str] | None = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        target_idx: Sequence[int] | None = None,
        n_samples: int = 1,
        **kwargs,
    ) -> pd.DataFrame:
        r"""Performs differential expression analysis.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`,
            defaults to the AnnData object used to initialize the model.
        groupby
            The key of the observations grouping to consider.
        group1
            Subset of groups, e.g. ["g1", "g2", "g3"], to which comparison shall be
            restricted, or all groups in `groupby` (default).
        group2
            If `None`, compare each group in `group1` to the union of the rest of
            the groups in `groupby`. If a group identifier, compare with respect to this
            group.
        idx1
            `idx1` and `idx2` can be used as an alternative to the AnnData keys.
            Custom identifier for `group1` that can be of three sorts:
            (1) a boolean mask, (2) indices, or (3) a string. If it is a string, then
            it will query indices that verifies conditions on adata.obs, as described
            in `pandas.DataFrame.query()`. If `idx1` is not `None`, this option
            overrides `group1` and `group2`.
        idx2
            Custom identifier for `group2` that has the same properties as `idx1`.
            By default, includes all cells not specified in `idx1`.
        mode:
            Method for differential expression. See
            https://docs.scvi-tools.org/en/0.14.1/user_guide/background/differential_expression.html
            for more details.
        delta
            Specific case of region inducing differential expression. In this case,
            we suppose that R\[-delta, delta] does not induce differential expression
            (change model default case).
        batch_size
            Mini-batch size for data loading into model. Defaults to
            scvi.settings.batch_size.
        all_stats
            Concatenate count statistics (e.g., mean expression group 1) to DE
            results.
        batch_correction
            Whether to correct for batch effects in DE inference.
        batchid1
            Subset of categories from `batch_key` registered in `setup_anndata`,
            e.g. ["batch1", "batch2", "batch3"], for `group1`. Only used if
            `batch_correction` is `True`, and by default all categories are used.
        batchid2
            Same as `batchid1` for `group2`. `batchid2` must either have null
            intersection with `batchid1`, or be exactly equal to `batchid1`. When the
            two sets are exactly equal, cells are compared by decoding on the same
            batch. When sets have null intersection, cells from `group1` and `group2`
            are decoded on each group in `group1` and `group2`, respectively.
        fdr_target
            Tag features as DE based on posterior expected false discovery rate.
        silent
            If `True`, disables the progress bar. Default: `False`.
        target_idx
            If not `None`, a boolean or integer identifier should be used for
            cells in the contrastive target group. Normalized expression values derived
            from both salient and background latent embeddings are used when
            {group1, group2} is a subset of the target group, otherwise background
            normalized expression values are used.
        kwargs: Keyword args for
            `scvi.model.base.DifferentialComputation.get_bayes_factors`.

        Returns
        -------
        Differential expression DataFrame.
        """
        adata = self._validate_anndata(adata)
        col_names = adata.var_names

        if target_idx is not None:
            target_idx = np.array(target_idx)
            if target_idx.dtype is np.dtype("bool"):
                assert len(target_idx) == adata.n_obs, (
                    "target_idx mask must be the same length as adata!"
                )
                target_idx = np.arange(adata.n_obs)[target_idx]
            model_fn = partial(
                self.get_specific_normalized_expression,
                return_numpy=True,
                n_samples=n_samples,
                batch_size=batch_size,
                expression_type=None,
                indices_to_return_salient=target_idx,
            )
        else:
            model_fn = partial(
                self.get_specific_normalized_expression,
                return_numpy=True,
                n_samples=n_samples,
                batch_size=batch_size,
                expression_type="salient",
                indices_to_return_salient=None,
            )

        result = _de_core(
            self.get_anndata_manager(adata, required=True),
            model_fn,
            representation_fn=None,
            groupby=groupby,
            group1=group1,
            group2=group2,
            idx1=idx1,
            idx2=idx2,
            all_stats=all_stats,
            all_stats_fn=scrna_raw_counts_properties,
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

    @staticmethod
    @torch.inference_mode()
    def _preprocess_normalized_expression(
        generative_outputs: dict[str, torch.Tensor],
        generative_output_key: str,
        gene_mask: list | slice,
        scaling: float,
    ) -> np.ndarray:
        output = generative_outputs[generative_output_key]
        output = output[..., gene_mask]
        output *= scaling
        output = output.cpu().numpy()
        return output

    @torch.inference_mode()
    def get_latent_library_size(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        give_mean: bool = True,
        batch_size: int | None = None,
    ) -> np.ndarray:
        r"""Returns the latent library size for each cell.

        This is denoted as :math:`\ell_n` in the scVI paper.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`,
            defaults to the AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        give_mean
            Return the mean or a sample from the posterior distribution.
        batch_size
            Minibatch size for data loading into model. Defaults to
            `scvi.settings.batch_size`.
        """
        self._check_if_trained(warn=False)

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        libraries = []
        for tensors in scdl:
            x = tensors[REGISTRY_KEYS.X_KEY]
            batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
            outputs = self.module._generic_inference(x=x, batch_index=batch_index)

            library = outputs["library"]
            if not give_mean:
                library = torch.exp(library)
            else:
                ql = (outputs["ql_m"], outputs["ql_v"])
                if ql is None:
                    raise RuntimeError(
                        "The module for this model does not compute the posterior"
                        "distribution for the library size. Set `give_mean` to False"
                        "to use the observed library size instead."
                    )
                library = torch.distributions.LogNormal(ql[0], ql[1]).mean
            libraries += [library.cpu()]
        return torch.cat(libraries).numpy()
