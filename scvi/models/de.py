import numpy as np
import inspect
from tqdm import tqdm
import pandas as pd
import torch
import logging
import warnings
from scvi.inference.posterior_utils import (
    pairs_sampler,
    describe_continuous_distrib,
    save_cluster_xlsx,
)
from typing import Union, List, Optional, Callable, Dict

logger = logging.getLogger(__name__)


class DifferentialExpression:
    def __init__(
        self,
        model_fn,  # method of model class (SCVI, not VAE)
        idx1: Union[List[bool], np.ndarray],
        idx2: Union[List[bool], np.ndarray],
        mode: Optional[str] = "vanilla",
        batchid1: Optional[Union[List[int], np.ndarray]] = None,
        batchid2: Optional[Union[List[int], np.ndarray]] = None,
        use_observed_batches: Optional[bool] = False,
        n_samples: int = 5000,
        use_permutation: bool = False,
        M_permutation: int = 10000,
        all_stats: bool = True,
        change_fn: Optional[Union[str, Callable]] = None,
        m1_domain_fn: Optional[Callable] = None,
        delta: Optional[float] = 0.5,
        cred_interval_lvls: Optional[Union[List[float], np.ndarray]] = None,
        **kwargs,
    ):
        r"""Unified method for differential expression inference.

        This function is an extension of the `get_bayes_factors` method
        providing additional genes information to the user

        Two modes coexist:

        - the "vanilla" mode follows protocol described in [Lopez18]_
        In this case, we perform hypothesis testing based on the hypotheses

        .. math::
            M_1: h_1 > h_2 ~\text{and}~ M_2: h_1 \leq h_2

        DE can then be based on the study of the Bayes factors

        .. math::
            \log p(M_1 | x_1, x_2) / p(M_2 | x_1, x_2)

        - the "change" mode (described in [Boyeau19]_)
        consists in estimating an effect size random variable (e.g., log fold-change) and
        performing Bayesian hypothesis testing on this variable.
        The `change_fn` function computes the effect size variable r based two inputs
        corresponding to the normalized means in both populations.

        Hypotheses:

        .. math::
            M_1: r \in R_1 ~\text{(effect size r in region inducing differential expression)}

        .. math::
            M_2: r  \notin R_1 ~\text{(no differential expression)}

        To characterize the region :math:`R_1`, which induces DE, the user has two choices.

        1. A common case is when the region :math:`[-\delta, \delta]` does not induce differential
        expression.
        If the user specifies a threshold delta,
        we suppose that :math:`R_1 = \mathbb{R} \setminus [-\delta, \delta]`

        2. specify an specific indicator function

        .. math::
            f: \mathbb{R} \mapsto \{0, 1\} ~\text{s.t.}~ r \in R_1 ~\text{iff.}~ f(r) = 1

        Decision-making can then be based on the estimates of

        .. math::
            p(M_1 \mid x_1, x_2)

        Both modes require to sample the normalized means posteriors.
        To that purpose, we sample the Posterior in the following way:

        1. The posterior is sampled n_samples times for each subpopulation

        2. For computation efficiency (posterior sampling is quite expensive), instead of
            comparing the obtained samples element-wise, we can permute posterior samples.
            Remember that computing the Bayes Factor requires sampling
            :math:`q(z_A \mid x_A)` and :math:`q(z_B \mid x_B)`

        Currently, the code covers several batch handling configurations:

        1. If ``use_observed_batches=True``, then batch are considered as observations
        and cells' normalized means are conditioned on real batch observations

        2. If case (cell group 1) and control (cell group 2) are conditioned on the same
        batch ids.
        Examples:

            >>> set(batchid1) = set(batchid2)

        or

            >>> batchid1 = batchid2 = None

        3. If case and control are conditioned on different batch ids that do not intersect
        i.e.,

            >>> set(batchid1) != set(batchid2)

        and

            >>> len(set(batchid1).intersection(set(batchid2))) == 0

        This function does not cover other cases yet and will warn users in such cases.

        Parameters
        ----------
        mode
            one of ["vanilla", "change"]
        idx1
            bool array masking subpopulation cells 1. Should be True where cell is
            from associated population
        idx2
            bool array masking subpopulation cells 2. Should be True where cell is
            from associated population
        batchid1
            List of batch ids for which you want to perform DE Analysis for
            subpopulation 1. By default, all ids are taken into account
        batchid2
            List of batch ids for which you want to perform DE Analysis for
            subpopulation 2. By default, all ids are taken into account
        use_observed_batches
            Whether normalized means are conditioned on observed
            batches
        n_samples
            Number of posterior samples
        use_permutation
            Activates step 2 described above.
            Simply formulated, pairs obtained from posterior sampling (when calling
            `sample_scale_from_batch`) will be randomly permuted so that the number of
            pairs used to compute Bayes Factors becomes M_permutation.
        M_permutation
            Number of times we will "mix" posterior samples in step 2.
            Only makes sense when use_permutation=True
        change_fn
            function computing effect size based on both normalized means
        m1_domain_fn
            custom indicator function of effect size regions
            inducing differential expression
        delta
            specific case of region inducing differential expression.
            In this case, we suppose that R \setminus [-\delta, \delta] does not induce differential expression
            (LFC case)
        cred_interval_lvls
            List of credible interval levels to compute for the posterior
            LFC distribution
        all_stats
            whether additional metrics should be provided
        **kwargs
            Other keywords arguments for `get_sample_scale`


        Returns
        -------
        diff_exp_results
            The most important columns are:

            - ``proba_de`` (probability of being differentially expressed in change mode)
            - ``bayes_factor`` (bayes factors in the vanilla mode)
            - ``scale1`` and ``scale2`` (means of the scales in population 1 and 2)
            - When using the change mode, the mean, median, std of the posterior LFC
        """

        all_info = self.get_bayes_factors(
            idx1=idx1,
            idx2=idx2,
            mode=mode,
            batchid1=batchid1,
            batchid2=batchid2,
            use_observed_batches=use_observed_batches,
            n_samples=n_samples,
            use_permutation=use_permutation,
            M_permutation=M_permutation,
            change_fn=change_fn,
            m1_domain_fn=m1_domain_fn,
            delta=delta,
            cred_interval_lvls=cred_interval_lvls,
            **kwargs,
        )
        gene_names = self.gene_dataset.gene_names
        if all_stats is True:
            (
                mean1,
                mean2,
                nonz1,
                nonz2,
                norm_mean1,
                norm_mean2,
            ) = self.gene_dataset.raw_counts_properties(idx1, idx2)
            genes_properties_dict = dict(
                raw_mean1=mean1,
                raw_mean2=mean2,
                non_zeros_proportion1=nonz1,
                non_zeros_proportion2=nonz2,
                raw_normalized_mean1=norm_mean1,
                raw_normalized_mean2=norm_mean2,
            )
            all_info = {**all_info, **genes_properties_dict}

        res = pd.DataFrame(all_info, index=gene_names)
        sort_key = "proba_de" if mode == "change" else "bayes_factor"
        res = res.sort_values(by=sort_key, ascending=False)
        return res

    def get_bayes_factors(
        self,
        model_fn,
        idx1: Union[List[bool], np.ndarray],
        idx2: Union[List[bool], np.ndarray],
        mode: Optional[str] = "vanilla",
        batchid1: Optional[Union[List[int], np.ndarray]] = None,
        batchid2: Optional[Union[List[int], np.ndarray]] = None,
        use_observed_batches: Optional[bool] = False,
        n_samples: int = 5000,
        use_permutation: bool = False,
        M_permutation: int = 10000,
        change_fn: Optional[Union[str, Callable]] = None,
        m1_domain_fn: Optional[Callable] = None,
        delta: Optional[float] = 0.5,
        cred_interval_lvls: Optional[Union[List[float], np.ndarray]] = None,
        **kwargs,
    ) -> Dict[str, np.ndarray]:
        r"""A unified method for differential expression inference.

        Two modes coexist:

        - the "vanilla" mode follows protocol described in [Lopez18]_
        In this case, we perform hypothesis testing based on the hypotheses

        .. math::
            M_1: h_1 > h_2 ~\text{and}~ M_2: h_1 \leq h_2

        DE can then be based on the study of the Bayes factors

        .. math::
            \log p(M_1 | x_1, x_2) / p(M_2 | x_1, x_2)

        - the "change" mode (described in [Boyeau19]_)
        consists in estimating an effect size random variable (e.g., log fold-change) and
        performing Bayesian hypothesis testing on this variable.
        The `change_fn` function computes the effect size variable r based two inputs
        corresponding to the normalized means in both populations.

        Hypotheses:

        .. math::
            M_1: r \in R_1 ~\text{(effect size r in region inducing differential expression)}

        .. math::
            M_2: r  \notin R_1 ~\text{(no differential expression)}

        To characterize the region :math:`R_1`, which induces DE, the user has two choices.

        1. A common case is when the region :math:`[-\delta, \delta]` does not induce differential
        expression.
        If the user specifies a threshold delta,
        we suppose that :math:`R_1 = \mathbb{R} \setminus [-\delta, \delta]`

        2. specify an specific indicator function

        .. math::
            f: \mathbb{R} \mapsto \{0, 1\} ~\text{s.t.}~ r \in R_1 ~\text{iff.}~ f(r) = 1

        Decision-making can then be based on the estimates of

        .. math::
            p(M_1 \mid x_1, x_2)

        Both modes require to sample the normalized means posteriors.
        To that purpose, we sample the Posterior in the following way:

        1. The posterior is sampled n_samples times for each subpopulation

        2. For computation efficiency (posterior sampling is quite expensive), instead of
            comparing the obtained samples element-wise, we can permute posterior samples.
            Remember that computing the Bayes Factor requires sampling
            :math:`q(z_A \mid x_A)` and :math:`q(z_B \mid x_B)`

        Currently, the code covers several batch handling configurations:

        1. If ``use_observed_batches=True``, then batch are considered as observations
        and cells' normalized means are conditioned on real batch observations

        2. If case (cell group 1) and control (cell group 2) are conditioned on the same
        batch ids.
        Examples:

            >>> set(batchid1) = set(batchid2)

        or

            >>> batchid1 = batchid2 = None


        3. If case and control are conditioned on different batch ids that do not intersect
        i.e.,
            >>> set(batchid1) != set(batchid2)

        and

            >>> len(set(batchid1).intersection(set(batchid2))) == 0


        This function does not cover other cases yet and will warn users in such cases.

        Parameters
        ----------
        mode
            one of ["vanilla", "change"]
        idx1
            bool array masking subpopulation cells 1. Should be True where cell is
            from associated population
        idx2
            bool array masking subpopulation cells 2. Should be True where cell is
            from associated population
        batchid1
            List of batch ids for which you want to perform DE Analysis for
            subpopulation 1. By default, all ids are taken into account
        batchid2
            List of batch ids for which you want to perform DE Analysis for
            subpopulation 2. By default, all ids are taken into account
        use_observed_batches
            Whether normalized means are conditioned on observed
            batches
        n_samples
            Number of posterior samples
        use_permutation
            Activates step 2 described above.
            Simply formulated, pairs obtained from posterior sampling (when calling
            `sample_scale_from_batch`) will be randomly permuted so that the number of
            pairs used to compute Bayes Factors becomes M_permutation.
        M_permutation
            Number of times we will "mix" posterior samples in step 2.
            Only makes sense when use_permutation=True
        change_fn
            function computing effect size based on both normalized means
        m1_domain_fn
            custom indicator function of effect size regions
            inducing differential expression
        delta
            specific case of region inducing differential expression.
            In this case, we suppose that :math:`R \setminus [-\delta, \delta]` does not induce differential expression
            (LFC case)
        cred_interval_lvls
            List of credible interval levels to compute for the posterior
            LFC distribution
        **kwargs:
            Other keywords arguments for `get_sample_scale()`


        Returns
        -------
        Differential expression properties

        """
        if not np.array_equal(self.indices, np.arange(len(self.gene_dataset))):
            logger.warning(
                "Differential expression requires a Posterior object created with all indices."
            )

        eps = 1e-8  # used for numerical stability
        # Normalized means sampling for both populations
        scales_batches_1 = self.scale_sampler(
            model_fn=model_fn,
            selection=idx1,
            batchid=batchid1,
            use_observed_batches=use_observed_batches,
            n_samples=n_samples,
            **kwargs,
        )
        scales_batches_2 = self.scale_sampler(
            model_fn=model_fn,
            selection=idx2,
            batchid=batchid2,
            use_observed_batches=use_observed_batches,
            n_samples=n_samples,
            **kwargs,
        )

        px_scale_mean1 = scales_batches_1["scale"].mean(axis=0)
        px_scale_mean2 = scales_batches_2["scale"].mean(axis=0)

        # Sampling pairs
        # The objective of code section below is to ensure than the samples of normalized
        # means we consider are conditioned on the same batch id
        batchid1_vals = np.unique(scales_batches_1["batch"])
        batchid2_vals = np.unique(scales_batches_2["batch"])

        create_pairs_from_same_batches = (
            set(batchid1_vals) == set(batchid2_vals)
        ) and not use_observed_batches
        if create_pairs_from_same_batches:
            # First case: same batch normalization in two groups
            logger.debug("Same batches in both cell groups")
            n_batches = len(set(batchid1_vals))
            n_samples_per_batch = (
                M_permutation // n_batches if M_permutation is not None else None
            )
            scales_1 = []
            scales_2 = []
            for batch_val in set(batchid1_vals):
                # Select scale samples that originate from the same batch id
                scales_1_batch = scales_batches_1["scale"][
                    scales_batches_1["batch"] == batch_val
                ]
                scales_2_batch = scales_batches_2["scale"][
                    scales_batches_2["batch"] == batch_val
                ]

                # Create more pairs
                scales_1_local, scales_2_local = pairs_sampler(
                    scales_1_batch,
                    scales_2_batch,
                    use_permutation=use_permutation,
                    M_permutation=n_samples_per_batch,
                )
                scales_1.append(scales_1_local)
                scales_2.append(scales_2_local)
            scales_1 = np.concatenate(scales_1, axis=0)
            scales_2 = np.concatenate(scales_2, axis=0)
        else:
            logger.debug("Ignoring batch conditionings to compare means")
            if len(set(batchid1_vals).intersection(set(batchid2_vals))) >= 1:
                warnings.warn(
                    "Batchids of cells groups 1 and 2 are different but have an non-null "
                    "intersection. Specific handling of such situations is not implemented "
                    "yet and batch correction is not trustworthy."
                )
            scales_1, scales_2 = pairs_sampler(
                scales_batches_1["scale"],
                scales_batches_2["scale"],
                use_permutation=use_permutation,
                M_permutation=M_permutation,
            )

        # Core of function: hypotheses testing based on the posterior samples we obtained above
        if mode == "vanilla":
            logger.debug("Differential expression using vanilla mode")
            proba_m1 = np.mean(scales_1 > scales_2, 0)
            proba_m2 = 1.0 - proba_m1
            res = dict(
                proba_m1=proba_m1,
                proba_m2=proba_m2,
                bayes_factor=np.log(proba_m1 + eps) - np.log(proba_m2 + eps),
                scale1=px_scale_mean1,
                scale2=px_scale_mean2,
            )

        elif mode == "change":
            logger.debug("Differential expression using change mode")

            # step 1: Construct the change function
            def lfc(x, y):
                return np.log2(x) - np.log2(y)

            if change_fn == "log-fold" or change_fn is None:
                change_fn = lfc
            elif not callable(change_fn):
                raise ValueError("'change_fn' attribute not understood")

            # step2: Construct the DE area function
            if m1_domain_fn is None:
                delta = delta if delta is not None else 0.5

                def m1_domain_fn(samples):
                    return np.abs(samples) >= delta

            change_fn_specs = inspect.getfullargspec(change_fn)
            domain_fn_specs = inspect.getfullargspec(m1_domain_fn)
            assert (len(change_fn_specs.args) == 2) & (
                len(domain_fn_specs.args) == 1
            ), "change_fn should take exactly two parameters as inputs; m1_domain_fn one parameter."
            try:
                change_distribution = change_fn(scales_1, scales_2)
                is_de = m1_domain_fn(change_distribution)
            except TypeError:
                raise TypeError(
                    "change_fn or m1_domain_fn have has wrong properties."
                    "Please ensure that these functions have the right signatures and"
                    "outputs and that they can process numpy arrays"
                )
            proba_m1 = np.mean(is_de, 0)
            change_distribution_props = describe_continuous_distrib(
                samples=change_distribution,
                credible_intervals_levels=cred_interval_lvls,
            )
            change_distribution_props = {
                "lfc_" + key: val for (key, val) in change_distribution_props.items()
            }

            res = dict(
                proba_de=proba_m1,
                proba_not_de=1.0 - proba_m1,
                bayes_factor=np.log(proba_m1 + eps) - np.log(1.0 - proba_m1 + eps),
                scale1=px_scale_mean1,
                scale2=px_scale_mean2,
                **change_distribution_props,
            )
        else:
            raise NotImplementedError("Mode {mode} not recognized".format(mode=mode))

        return res

    @torch.no_grad()
    def scale_sampler(
        self,
        model_fn,
        selection: Union[List[bool], np.ndarray],
        n_samples: Optional[int] = 5000,
        n_samples_per_cell: Optional[int] = None,
        batchid: Optional[Union[List[int], np.ndarray]] = None,
        use_observed_batches: Optional[bool] = False,
        give_mean: Optional[bool] = False,
        **kwargs,
    ) -> dict:
        """Samples the posterior scale using the variational posterior distribution.

        Parameters
        ----------
        n_samples
            Number of samples in total per batch (fill either `n_samples_total`
            or `n_samples_per_cell`)
        n_samples_per_cell
            Number of time we sample from each observation per batch
            (fill either `n_samples_total` or `n_samples_per_cell`)
        batchid
            Biological batch for which to sample from.
            Default (None) sample from all batches
        use_observed_batches
            Whether normalized means are conditioned on observed
            batches or if observed batches are to be used
        selection
            Mask or list of cell ids to select
        **kwargs
            Other keywords arguments for `get_sample_scale()`


        Returns
        -------
        type
            Dictionary containing:
            `scale`
            Posterior aggregated scale samples of shape (n_samples, n_genes)
            where n_samples correspond to either:
            - n_bio_batches * n_cells * n_samples_per_cell
            or
            - n_samples_total
            `batch`
            associated batch ids
        """
        # Get overall number of desired samples and desired batches
        if batchid is None and not use_observed_batches:
            batchid = np.arange(self.gene_dataset.n_batches)
        if use_observed_batches:
            assert batchid is None, "Unconsistent batch policy"
            batchid = [None]
        if n_samples is None and n_samples_per_cell is None:
            n_samples = 5000
        elif n_samples_per_cell is not None and n_samples is None:
            n_samples = n_samples_per_cell * len(selection)
        if (n_samples_per_cell is not None) and (n_samples is not None):
            warnings.warn(
                "n_samples and n_samples_per_cell were provided. Ignoring n_samples_per_cell"
            )
        n_samples = int(n_samples / len(batchid))
        if n_samples == 0:
            warnings.warn(
                "very small sample size, please consider increasing `n_samples`"
            )
            n_samples = 2

        # Selection of desired cells for sampling
        if selection is None:
            raise ValueError("selections should be a list of cell subsets indices")
        selection = np.array(selection)
        if selection.dtype is np.dtype("bool"):
            selection = np.asarray(np.where(selection)[0].ravel())
        old_loader = self.data_loader

        # Sampling loop
        px_scales = []
        batch_ids = []
        for batch_idx in batchid:
            idx = np.random.choice(
                np.arange(len(self.gene_dataset))[selection], n_samples
            )
            self.update_sampler_indices(idx=idx)
            px_scales.append(model_fn(transform_batch=batch_idx, **kwargs))
            batch_idx = batch_idx if batch_idx is not None else np.nan
            batch_ids.append(np.ones((px_scales[-1].shape[0])) * batch_idx)
        px_scales = np.concatenate(px_scales)
        batch_ids = np.concatenate(batch_ids).reshape(-1)
        assert (
            px_scales.shape[0] == batch_ids.shape[0]
        ), "sampled scales and batches have inconsistent shapes"
        self.data_loader = old_loader
        if give_mean:
            px_scales = px_scales.mean(0)
        return dict(scale=px_scales, batch=batch_ids)

    @torch.no_grad()
    def one_vs_all_degenes(
        self,
        subset: Optional[Union[List[bool], np.ndarray]] = None,
        cell_labels: Optional[Union[List, np.ndarray]] = None,
        use_observed_batches: bool = False,
        min_cells: int = 10,
        n_samples: int = 5000,
        use_permutation: bool = False,
        M_permutation: int = 10000,
        output_file: bool = False,
        mode: Optional[str] = "vanilla",
        change_fn: Optional[Union[str, Callable]] = None,
        m1_domain_fn: Optional[Callable] = None,
        delta: Optional[float] = 0.5,
        cred_interval_lvls: Optional[Union[List[float], np.ndarray]] = None,
        save_dir: str = "./",
        filename="one2all",
        **kwargs,
    ) -> tuple:
        r"""Performs one population vs all others Differential Expression Analysis

        It takes labels or cell types to characterize the different populations.

        Parameters
        ----------
        subset
            None Or bool array masking subset of cells you are interested in
            (True when you want to select cell). In that case, it should have same length than `gene_dataset`
        cell_labels
            optional: Labels of cells
        min_cells
            Ceil number of cells used to compute Bayes Factors
        n_samples
            Number of times the posterior will be sampled for each pop
        use_permutation
            Activates pair random permutations.
            Simply formulated, pairs obtained from posterior sampling (when calling
            `sample_scale_from_batch`) will be randomly permuted so that the number of
            pairs used to compute Bayes Factors becomes M_permutation.
        M_permutation
            Number of times we will "mix" posterior samples in step 2.
            Only makes sense when use_permutation=True
        use_observed_batches
            see `differential_expression_score`
        M_permutation
            see `differential_expression_score`
        mode
            see `differential_expression_score`
        change_fn
            see `differential_expression_score`
        m1_domain_fn
            see `differential_expression_score`
        delta
            see `differential_expression_score
        cred_interval_lvls
            List of credible interval levels to compute for the posterior
            LFC distribution
        output_file
            Bool: save file?
        save_dir
            param filename:`
        **kwargs
            Other keywords arguments for `get_sample_scale`

        Returns
        -------
        type
            Tuple (de_res, de_cluster) (i) de_res is a list of length nb_clusters
            (based on provided labels or on hardcoded cell types) (ii) de_res[i] contains Bayes Factors
            for population number i vs all the rest (iii) de_cluster returns the associated names of clusters.
            Are contained in this results only clusters for which we have at least `min_cells`
            elements to compute predicted Bayes Factors
        """
        if cell_labels is not None:
            if len(cell_labels) != len(self.gene_dataset):
                raise ValueError(
                    " the length of cell_labels have to be "
                    "the same as the number of cells"
                )
        # Input cell_labels take precedence over cell type label annotation in dataset
        if cell_labels is not None:
            cluster_id = np.unique(cell_labels[cell_labels >= 0])
            # Can make cell_labels < 0 to filter out cells when computing DE
        else:
            cell_labels = self.gene_dataset.labels.ravel()
            cluster_id = np.unique(cell_labels)

        de_res = []
        de_cluster = []
        for i, x in enumerate(tqdm(cluster_id)):
            if subset is None:
                idx1 = cell_labels == i
                idx2 = cell_labels != i
            else:
                idx1 = (cell_labels == i) * subset
                idx2 = (cell_labels != i) * subset
            if np.sum(idx1) > min_cells and np.sum(idx2) > min_cells:
                de_cluster.append(x)
                res = self.differential_expression_score(
                    idx1=idx1,
                    idx2=idx2,
                    mode=mode,
                    change_fn=change_fn,
                    m1_domain_fn=m1_domain_fn,
                    delta=delta,
                    use_observed_batches=use_observed_batches,
                    M_permutation=M_permutation,
                    n_samples=n_samples,
                    use_permutation=use_permutation,
                    cred_interval_lvls=cred_interval_lvls,
                    **kwargs,
                )
                res["clusters"] = np.repeat(x, len(res.index))
                de_res.append(res)
        if output_file:  # store as an excel spreadsheet
            save_cluster_xlsx(
                filepath=save_dir + "differential_expression.%s.xlsx" % filename,
                cluster_names=de_cluster,
                de_results=de_res,
            )
        return de_res, de_cluster

    def within_cluster_degenes(
        self,
        states: Union[List[bool], np.ndarray],
        cell_labels: Optional[Union[List, np.ndarray]] = None,
        min_cells: int = 10,
        batch1: Optional[Union[List[int], np.ndarray]] = None,
        batch2: Optional[Union[List[int], np.ndarray]] = None,
        use_observed_batches: bool = False,
        subset: Optional[Union[List[bool], np.ndarray]] = None,
        n_samples: int = 5000,
        use_permutation: bool = False,
        M_permutation: int = 10000,
        mode: Optional[str] = "vanilla",
        change_fn: Optional[Union[str, Callable]] = None,
        m1_domain_fn: Optional[Callable] = None,
        delta: Optional[float] = 0.5,
        cred_interval_lvls: Optional[Union[List[float], np.ndarray]] = None,
        output_file: bool = False,
        save_dir: str = "./",
        filename: str = "within_cluster",
        **kwargs,
    ) -> tuple:
        """Performs Differential Expression within clusters for different cell states

        Parameters
        ----------
        cell_labels
            optional: Labels of cells
        min_cells
            Ceil number of cells used to compute Bayes Factors
        states
            States of the cells.
        batch1
            List of batch ids for which you want to perform DE Analysis for
            subpopulation 1. By default, all ids are taken into account
        batch2
            List of batch ids for which you want to perform DE Analysis for
            subpopulation 2. By default, all ids are taken into account
        subset
            MASK: Subset of cells you are interested in.
        n_samples
            Number of times the posterior will be sampled for each pop
        use_permutation
            Activates pair random permutations.
            Simply formulated, pairs obtained from posterior sampling (when calling
            `sample_scale_from_batch`) will be randomly permuted so that the number of
            pairs used to compute Bayes Factors becomes M_permutation.
        M_permutation
            Number of times we will "mix" posterior samples in step 2.
            Only makes sense when use_permutation=True
        output_file
            Bool: save file?
        save_dir
            param filename:
        use_observed_batches
            see `differential_expression_score`
        M_permutation
            see `differential_expression_score`
        mode
            see `differential_expression_score`
        change_fn
            see `differential_expression_score`
        m1_domain_fn
            see `differential_expression_score`
        delta
            see `differential_expression_score`
        cred_interval_lvls
            See `differential_expression_score`
        **kwargs
            Other keywords arguments for `get_sample_scale()`

        Returns
        -------
        type
            Tuple (de_res, de_cluster) (i) de_res is a list of length nb_clusters
            (based on provided labels or on hardcoded cell types) (ii) de_res[i] contains Bayes Factors
            for population number i vs all the rest (iii) de_cluster returns the associated names of clusters.
            Are contained in this results only clusters for which we have at least `min_cells`
            elements to compute predicted Bayes Factors
        """
        if len(self.gene_dataset) != len(states):
            raise ValueError(
                " the length of states have to be the same as the number of cells"
            )
        if cell_labels is not None:
            if len(cell_labels) != len(self.gene_dataset):
                raise ValueError(
                    " the length of cell_labels have to be "
                    "the same as the number of cells"
                )
        if (cell_labels is None) and not hasattr(self.gene_dataset, "cell_types"):
            raise ValueError(
                "If gene_dataset is not annotated with labels and cell types,"
                " then must provide cell_labels"
            )
        # Input cell_labels take precedence over cell type label annotation in dataset
        elif cell_labels is not None:
            cluster_id = np.unique(cell_labels[cell_labels >= 0])
            # Can make cell_labels < 0 to filter out cells when computing DE
        else:
            cell_labels = self.gene_dataset.labels.ravel()
            cluster_id = np.unique(cell_labels)
        de_res = []
        de_cluster = []
        oppo_states = ~states
        for i, x in enumerate(tqdm(cluster_id)):
            if subset is None:
                idx1 = (cell_labels == i) * states
                idx2 = (cell_labels == i) * oppo_states
            else:
                idx1 = (cell_labels == i) * subset * states
                idx2 = (cell_labels == i) * subset * oppo_states
            if np.sum(idx1) > min_cells and np.sum(idx2) > min_cells:
                de_cluster.append(x)
                res = self.differential_expression_score(
                    idx1=idx1.astype(bool),
                    idx2=idx2.astype(bool),
                    batchid1=batch1,
                    batchid2=batch2,
                    use_observed_batches=use_observed_batches,
                    M_permutation=M_permutation,
                    n_samples=n_samples,
                    use_permutation=use_permutation,
                    mode=mode,
                    change_fn=change_fn,
                    m1_domain_fn=m1_domain_fn,
                    delta=delta,
                    cred_interval_lvls=cred_interval_lvls,
                    **kwargs,
                )
                res["clusters"] = np.repeat(x, len(res.index))
                de_res.append(res)
        if output_file:  # store as an excel spreadsheet
            save_cluster_xlsx(
                filepath=save_dir + "differential_expression.%s.xlsx" % filename,
                cluster_names=de_cluster,
                de_results=de_res,
            )
        return de_res, de_cluster
