from functools import partial
from typing import Literal

import joblib
import pandas as pd
from anndata import AnnData

from scvi.model._utils import scrna_raw_counts_properties
from scvi.model.base import (
    RNASeqMixin,
)
from scvi.model.base._de_core import _de_core
from scvi.utils import de_dsp

from .differential_expression import _niche_de_core


class NicheRNASeqMixin(RNASeqMixin):
    @de_dsp.dedent
    def differential_expression(
        self,
        adata: AnnData | None = None,
        groupby: str | None = None,
        group1: list[str] | None = None,
        group2: str | None = None,
        idx1: list[int] | list[bool] | str | None = None,
        idx2: list[int] | list[bool] | str | None = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float | list[float] = 0.25,
        batch_size: int | None = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: list[str] | None = None,
        batchid2: list[str] | None = None,
        fdr_target: float | list[float] = 0.05,
        silent: bool = False,
        weights: Literal["uniform", "importance"] | None = "uniform",
        filter_outlier_cells: bool = False,
        importance_weighting_kwargs: dict | None = None,
        ###### NicheSCVI specific ######
        radius: int | None = 50,
        k_nn: int | None = None,
        niche_mode: bool = True,
        n_restarts_optimizer_gpc: int = 10,
        path_to_save: str | None = None,
        **kwargs,
    ) -> pd.DataFrame:
        r"""A unified method for differential expression analysis.

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
        weights
            Weights to use for sampling. If `None`, defaults to `"uniform"`.
        filter_outlier_cells
            Whether to filter outlier cells with
            :meth:`~scvi.model.base.DifferentialComputation.filter_outlier_cells`.
        importance_weighting_kwargs
            Keyword arguments passed into
            :meth:`~scvi.model.base.RNASeqMixin._get_importance_weights`.
        radius
            Radius for NicheSCVI DE.
        k_nn
            Number of nearest neighbors for NicheSCVI DE.
        niche_mode
            Whether to use NicheSCVI DE or SCVI DE.
        **kwargs
            Keyword args for :meth:`scvi.model.base.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential expression DataFrame.
        """
        adata = self._validate_anndata(adata)
        col_names = adata.var_names
        importance_weighting_kwargs = importance_weighting_kwargs or {}
        model_fn = partial(
            self.get_normalized_expression,
            return_numpy=True,
            n_samples=1,
            batch_size=batch_size,
            weights=weights,
            **importance_weighting_kwargs,
        )
        representation_fn = self.get_latent_representation if filter_outlier_cells else None

        if niche_mode:
            result = _niche_de_core(
                self.get_anndata_manager(adata, required=True),
                model_fn,
                representation_fn,
                groupby,
                group1,
                group2,
                idx1,
                idx2,
                all_stats,
                scrna_raw_counts_properties,
                col_names,
                mode,
                batchid1,
                batchid2,
                delta,
                batch_correction,
                fdr_target,
                silent,
                radius=radius,
                k_nn=k_nn,
                n_restarts_optimizer_gpc=n_restarts_optimizer_gpc,
                **kwargs,
            )

        else:
            result = _de_core(
                self.get_anndata_manager(adata, required=True),
                model_fn,
                representation_fn,
                groupby,
                group1,
                group2,
                idx1,
                idx2,
                all_stats,
                scrna_raw_counts_properties,
                col_names,
                mode,
                batchid1,
                batchid2,
                delta,
                batch_correction,
                fdr_target,
                silent,
                **kwargs,
            )

        if path_to_save is not None:
            joblib.dump(result, path_to_save)

        return result
