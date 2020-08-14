import numpy as np
import pandas as pd
import torch

from typing import List, Optional, Union


def pairs_sampler(
    arr1: Union[List[float], np.ndarray, torch.Tensor],
    arr2: Union[List[float], np.ndarray, torch.Tensor],
    use_permutation: bool = True,
    M_permutation: int = None,
    sanity_check_perm: bool = False,
    weights1: Union[List[float], np.ndarray, torch.Tensor] = None,
    weights2: Union[List[float], np.ndarray, torch.Tensor] = None,
) -> tuple:
    """In a context where we want to estimate a double sum, virtually increases the number
    of samples by considering more pairs so as to better estimate the double summation operation

    Parameters
    ----------
    arr1
        samples from population 1
    arr2
        samples from population 2
    use_permutation
        Whether to mix samples from both populations
    M_permutation
        param sanity_check_perm: If True, resulting mixed arrays arr1 and arr2 are mixed together
        In most cases, this parameter should remain False
    weights1
        probabilities associated to array 1 for random sampling
    weights2
        probabilities associated to array 2 for random sampling

    Returns
    -------
    type
        new_arr1, new_arr2

    """
    if use_permutation is True:
        # prepare the pairs for sampling
        n_arr1 = arr1.shape[0]
        n_arr2 = arr2.shape[0]
        if not sanity_check_perm:
            # case1: no permutation, sample from A and then from B
            u, v = (
                np.random.choice(n_arr1, size=M_permutation, p=weights1),
                np.random.choice(n_arr2, size=M_permutation, p=weights2),
            )
            first_set = arr1[u]
            second_set = arr2[v]
        else:
            # case2: permutation, sample from A+B twice (sanity check)
            u, v = (
                np.random.choice(n_arr1 + n_arr2, size=M_permutation),
                np.random.choice(n_arr1 + n_arr2, size=M_permutation),
            )
            concat_arr = np.concatenate((arr1, arr2))
            first_set = concat_arr[u]
            second_set = concat_arr[v]
    else:
        first_set = arr1
        second_set = arr2
    return first_set, second_set


def credible_intervals(
    ary: np.ndarray, confidence_level: Union[float, List[float], np.ndarray] = 0.94
) -> np.ndarray:
    """Calculate highest posterior density (HPD) of array for given credible_interval.

    Taken from the arviz package
    The HPD is the minimum width Bayesian credible interval (BCI). This implementation works only
    for unimodal distributions.

    Parameters
    ----------
    ary
        posterior samples
    confidence_level
        confidence level

    Returns
    -------
    type
        intervals minima, intervals maxima

    """
    if ary.ndim > 1:
        hpd = np.array(
            [
                credible_intervals(row, confidence_level=confidence_level)
                for row in ary.T
            ]
        )
        return hpd
    # Make a copy of trace
    ary = ary.copy()
    n = len(ary)
    ary = np.sort(ary)
    interval_idx_inc = int(np.floor(confidence_level * n))
    n_intervals = n - interval_idx_inc
    interval_width = ary[interval_idx_inc:] - ary[:n_intervals]

    if len(interval_width) == 0:
        raise ValueError(
            "Too few elements for interval calculation. "
            "Check that credible_interval meets condition 0 =< credible_interval < 1"
        )
    min_idx = np.argmin(interval_width)
    hdi_min = ary[min_idx]
    hdi_max = ary[min_idx + interval_idx_inc]
    return np.array([hdi_min, hdi_max])


def describe_continuous_distrib(
    samples: Union[np.ndarray, torch.Tensor],
    credible_intervals_levels: Optional[Union[List[float], np.ndarray]] = None,
) -> dict:
    """Computes properties of distribution based on its samples

    Parameters
    ----------
    samples
        samples of shape (n_samples, n_features)
    credible_intervals_levels
        Confidence in (0, 1)
        of credible intervals to be computed

    Returns
    -------
    type
        properties of distribution

    """
    dist_props = dict(
        mean=samples.mean(0),
        median=np.median(samples, 0),
        std=samples.std(0),
        min=samples.min(0),
        max=samples.max(0),
    )
    credible_intervals_levels = (
        [] if credible_intervals_levels is None else credible_intervals_levels
    )
    for confidence in credible_intervals_levels:
        intervals = credible_intervals(samples, confidence_level=confidence)
        interval_min, interval_max = intervals[:, 0], intervals[:, 1]
        conf_str = str(confidence)[:5]
        dist_props["confidence_interval_{}_min".format(conf_str)] = interval_min
        dist_props["confidence_interval_{}_max".format(conf_str)] = interval_max

    return dist_props


def save_cluster_xlsx(
    filepath: str, de_results: List[pd.DataFrame], cluster_names: List
):
    """Saves multi-clusters DE in an xlsx sheet

    Parameters
    ----------
    filepath
        xslx save path
    de_results
        list of pandas Dataframes for each cluster
    cluster_names
        list of cluster names
    """
    writer = pd.ExcelWriter(filepath, engine="xlsxwriter")
    for i, x in enumerate(cluster_names):
        de_results[i].to_excel(writer, sheet_name=str(x))
    writer.close()
