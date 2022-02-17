import logging
from functools import partial
from typing import Iterable, Optional, Sequence, Union

import numpy as np
import pandas as pd
import torch
from anndata import AnnData
from sklearn.covariance import EllipticEnvelope
from torch.distributions import Categorical, Normal

from scvi._compat import Literal
from scvi._utils import _doc_params
from scvi.model._utils import _get_batch_code_from_category, scrna_raw_counts_properties
from scvi.model.base._utils import _de_core
from scvi.utils._docstrings import doc_differential_expression

logger = logging.getLogger(__name__)
Number = Union[int, float]


class DEMixin:
    """
    DE module relying on Importance-sampling.

    Mixin for using
    importance-weighted DE content.
    This however requires some additional structure on the
    module's (e.g., VAE) methods and associate signatures
    """

    @_doc_params(
        doc_differential_expression=doc_differential_expression,
    )
    def differential_expression(
        self,
        adata: Optional[AnnData] = None,
        groupby: Optional[str] = None,
        group1: Optional[Iterable[str]] = None,
        group2: Optional[str] = None,
        idx1: Optional[Union[Sequence[int], Sequence[bool], str]] = None,
        idx2: Optional[Union[Sequence[int], Sequence[bool], str]] = None,
        mode: Literal["vanilla", "change"] = "change",
        delta: float = 0.25,
        batch_size: Optional[int] = None,
        all_stats: bool = True,
        batch_correction: bool = False,
        batchid1: Optional[Iterable[str]] = None,
        batchid2: Optional[Iterable[str]] = None,
        fdr_target: float = 0.05,
        silent: bool = False,
        pseudocounts: float = 0.0,
        fn_kwargs: Optional[dict] = None,
        importance_sampling: Optional[bool] = False,
        **kwargs,
    ) -> pd.DataFrame:
        r"""
        A unified method for differential expression analysis.

        Implements `"vanilla"` DE [Lopez18]_ and `"change"` mode DE [Boyeau19]_.
        When using the change method, uses either the plugin estimator
        or importance sampling for improved FDR control.

        Parameters
        ----------
        {doc_differential_expression}
        **kwargs
            Keyword args for :meth:`scvi.model.base.DifferentialComputation.get_bayes_factors`

        Returns
        -------
        Differential expression DataFrame.
        """
        adata = self._validate_anndata(adata)

        col_names = adata.var_names
        model_fn = partial(
            self.get_normalized_expression,
            return_numpy=True,
            n_samples=1,
            batch_size=batch_size,
        )
        result = _de_core(
            self.get_anndata_manager(adata, required=True),
            model_fn,
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
            pseudocounts=pseudocounts,
            use_permutation=True,
            **kwargs,
        )

        return result
