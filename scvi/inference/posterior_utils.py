import os
from typing import List, Optional, Union

import anndata
import numpy as np
import pandas as pd
import scipy
import torch
import torch.nn as nn
from matplotlib import pyplot as plt
from scipy.optimize import linear_sum_assignment
from scipy.stats import entropy, kde
from sklearn.neighbors import KNeighborsRegressor, NearestNeighbors

from scvi.dataset import AnnDatasetFromAnnData


def load_posterior(
    dir_path: str,
    model: nn.Module,
    use_cuda: Optional[Union[bool, str]] = "auto",
    **posterior_kwargs
):
    """Function to use in order to retrieve a posterior that was saved using the ``save_posterior`` method

    Because of pytorch model loading usage, this function needs a scVI model object initialized with exact same parameters
    that during training.
    Because saved posteriors correspond to already trained models, data is loaded sequentially using a ``SequentialSampler``.

    :param dir_path: directory containing the posterior properties to be retrieved.
    :param model: scVI initialized model.
    :param use_cuda: Specifies if the computations should be perfomed with a GPU.
      Default: ``True``
      If ``auto``, then cuda availability is inferred, with a preference to load on GPU.
      If ``False``, the model will be loaded on the CPU, even if it was trained using a GPU.
    :param posterior_kwargs: additional parameters to feed to the posterior constructor.


    Usage example:
    1. Save posterior
        >>> model = VAE(nb_genes, n_batches, n_hidden=128, n_latent=10)
        >>> trainer = UnsupervisedTrainer(vae, dataset, train_size=0.5, use_cuda=use_cuda)
        >>> trainer.train(n_epochs=200)
        >>> trainer.train_set.save_posterior("./my_run_train_posterior")

    2. Load posterior
        >>> model = VAE(nb_genes, n_batches, n_hidden=128, n_latent=10)
        >>> post = load_posterior("./my_run_train_posterior", model=model)
    """
    # Avoid circular imports
    from scvi.inference.total_inference import TotalPosterior
    from scvi.inference.jvae_trainer import JPosterior
    from scvi.inference.posterior import Posterior
    from scvi.inference.annotation import AnnotationPosterior

    post_type_path = os.path.join(dir_path, "posterior_type.txt")
    dataset_path = os.path.join(dir_path, "anndata_dataset.h5ad")
    model_path = os.path.join(dir_path, "model_params.pt")
    indices_path = os.path.join(dir_path, "indices.npy")
    data_loader_kwargs_path = os.path.join(dir_path, "data_loader_kwargs.h5")

    # Infering posterior type
    with open(post_type_path, "r") as post_file:
        post_class_str = post_file.readline()
    str_to_classes = dict(
        TotalPosterior=TotalPosterior,
        JPosterior=JPosterior,
        Posterior=Posterior,
        AnnotationPosterior=AnnotationPosterior,
    )
    if post_class_str not in str_to_classes:
        raise ValueError(
            "Posterior type {} not eligible for loading".format(post_class_str)
        )
    post_class = str_to_classes[post_class_str]

    # Loading dataset and associated measurements
    ad = anndata.read_h5ad(filename=dataset_path)
    key = "cell_measurements_col_mappings"
    if key in ad.uns:
        cell_measurements_col_mappings = ad.uns[key]
    else:
        cell_measurements_col_mappings = dict()
    dataset = AnnDatasetFromAnnData(
        ad=ad, cell_measurements_col_mappings=cell_measurements_col_mappings
    )

    # Loading scVI model
    if use_cuda == "auto":
        use_cuda = torch.cuda.is_available()
    use_cuda = use_cuda and torch.cuda.is_available()
    if use_cuda:
        model.load_state_dict(torch.load(model_path))
        model.cuda()
    else:
        device = torch.device("cpu")
        model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()

    # Loading data loader options and posterior
    indices = np.load(file=indices_path)
    data_loader_kwargs = pd.read_hdf(
        data_loader_kwargs_path, key="data_loader"
    ).to_dict()
    my_post = post_class(
        model=model,
        gene_dataset=dataset,
        shuffle=False,
        indices=indices,
        use_cuda=use_cuda,
        data_loader_kwargs=data_loader_kwargs,
        **posterior_kwargs
    )
    return my_post


def entropy_from_indices(indices):
    return entropy(np.array(np.unique(indices, return_counts=True)[1].astype(np.int32)))


def entropy_batch_mixing(
    latent_space, batches, n_neighbors=50, n_pools=50, n_samples_per_pool=100
):
    def entropy(hist_data):
        n_batches = len(np.unique(hist_data))
        if n_batches > 2:
            raise ValueError("Should be only two clusters for this metric")
        frequency = np.mean(hist_data == 1)
        if frequency == 0 or frequency == 1:
            return 0
        return -frequency * np.log(frequency) - (1 - frequency) * np.log(1 - frequency)

    n_neighbors = min(n_neighbors, len(latent_space) - 1)
    nne = NearestNeighbors(n_neighbors=1 + n_neighbors, n_jobs=8)
    nne.fit(latent_space)
    kmatrix = nne.kneighbors_graph(latent_space) - scipy.sparse.identity(
        latent_space.shape[0]
    )

    score = 0
    for t in range(n_pools):
        indices = np.random.choice(
            np.arange(latent_space.shape[0]), size=n_samples_per_pool
        )
        score += np.mean(
            [
                entropy(
                    batches[
                        kmatrix[indices].nonzero()[1][
                            kmatrix[indices].nonzero()[0] == i
                        ]
                    ]
                )
                for i in range(n_samples_per_pool)
            ]
        )
    return score / float(n_pools)


def plot_imputation(original, imputed, show_plot=True, title="Imputation"):
    y = imputed
    x = original

    ymax = 10
    mask = x < ymax
    x = x[mask]
    y = y[mask]

    mask = y < ymax
    x = x[mask]
    y = y[mask]

    l_minimum = np.minimum(x.shape[0], y.shape[0])

    x = x[:l_minimum]
    y = y[:l_minimum]

    data = np.vstack([x, y])

    plt.figure(figsize=(5, 5))

    axes = plt.gca()
    axes.set_xlim([0, ymax])
    axes.set_ylim([0, ymax])

    nbins = 50

    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    k = kde.gaussian_kde(data)
    xi, yi = np.mgrid[0 : ymax : nbins * 1j, 0 : ymax : nbins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    plt.title(title, fontsize=12)
    plt.ylabel("Imputed counts")
    plt.xlabel("Original counts")

    plt.pcolormesh(yi, xi, zi.reshape(xi.shape), cmap="Reds")

    a, _, _, _ = np.linalg.lstsq(y[:, np.newaxis], x, rcond=-1)
    linspace = np.linspace(0, ymax)
    plt.plot(linspace, a * linspace, color="black")

    plt.plot(linspace, linspace, color="black", linestyle=":")
    if show_plot:
        plt.show()
    plt.savefig(title + ".png")


def nn_overlap(X1, X2, k=100):
    """
    Compute the overlap between the k-nearest neighbor graph of X1 and X2 using Spearman correlation of the
    adjacency matrices.
    """
    assert len(X1) == len(X2)
    n_samples = len(X1)
    k = min(k, n_samples - 1)
    nne = NearestNeighbors(n_neighbors=k + 1)  # "n_jobs=8
    nne.fit(X1)
    kmatrix_1 = nne.kneighbors_graph(X1) - scipy.sparse.identity(n_samples)
    nne.fit(X2)
    kmatrix_2 = nne.kneighbors_graph(X2) - scipy.sparse.identity(n_samples)

    # 1 - spearman correlation from knn graphs
    spearman_correlation = scipy.stats.spearmanr(
        kmatrix_1.A.flatten(), kmatrix_2.A.flatten()
    )[0]
    # 2 - fold enrichment
    set_1 = set(np.where(kmatrix_1.A.flatten() == 1)[0])
    set_2 = set(np.where(kmatrix_2.A.flatten() == 1)[0])
    fold_enrichment = (
        len(set_1.intersection(set_2))
        * n_samples ** 2
        / (float(len(set_1)) * len(set_2))
    )
    return spearman_correlation, fold_enrichment


def unsupervised_clustering_accuracy(
    y: Union[np.ndarray, torch.Tensor], y_pred: Union[np.ndarray, torch.Tensor]
) -> tuple:
    """
    Unsupervised Clustering Accuracy
    """
    assert len(y_pred) == len(y)
    u = np.unique(np.concatenate((y, y_pred)))
    n_clusters = len(u)
    mapping = dict(zip(u, range(n_clusters)))
    reward_matrix = np.zeros((n_clusters, n_clusters), dtype=np.int64)
    for y_pred_, y_ in zip(y_pred, y):
        if y_ in mapping:
            reward_matrix[mapping[y_pred_], mapping[y_]] += 1
    cost_matrix = reward_matrix.max() - reward_matrix
    row_assign, col_assign = linear_sum_assignment(cost_matrix)

    # Construct optimal assignments matrix
    row_assign = row_assign.reshape((-1, 1))  # (n,) to (n, 1) reshape
    col_assign = col_assign.reshape((-1, 1))  # (n,) to (n, 1) reshape
    assignments = np.concatenate((row_assign, col_assign), axis=1)

    optimal_reward = reward_matrix[row_assign, col_assign].sum() * 1.0
    return optimal_reward / y_pred.size, assignments


def knn_purity(latent, label, n_neighbors=30):
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(latent)
    indices = nbrs.kneighbors(latent, return_distance=False)[:, 1:]
    neighbors_labels = np.vectorize(lambda i: label[i])(indices)

    # pre cell purity scores
    scores = ((neighbors_labels - label.reshape(-1, 1)) == 0).mean(axis=1)
    res = [
        np.mean(scores[label == i]) for i in np.unique(label)
    ]  # per cell-type purity

    return np.mean(res)


def proximity_imputation(real_latent1, normed_gene_exp_1, real_latent2, k=4):
    knn = KNeighborsRegressor(k, weights="distance")
    y = knn.fit(real_latent1, normed_gene_exp_1).predict(real_latent2)
    return y


def pairs_sampler(
    arr1: Union[List[float], np.ndarray, torch.Tensor],
    arr2: Union[List[float], np.ndarray, torch.Tensor],
    use_permutation: bool = True,
    M_permutation: int = None,
    sanity_check_perm: bool = False,
    weights1: Union[List[float], np.ndarray, torch.Tensor] = None,
    weights2: Union[List[float], np.ndarray, torch.Tensor] = None,
) -> tuple:
    """
    In a context where we want to estimate a double sum, virtually increases the number
    of samples by considering more pairs so as to better estimate the double summation operation

    :param arr1: samples from population 1
    :param arr2: samples from population 2
    :param use_permutation: Whether to mix samples from both populations
    :param M_permutation:
    :param sanity_check_perm: If True, resulting mixed arrays arr1 and arr2 are mixed together
    In most cases, this parameter should remain False
    :param weights1: probabilities associated to array 1 for random sampling
    :param weights2: probabilities associated to array 2 for random sampling
    :return: new_arr1, new_arr2
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
    """
    Taken from the arviz package
    Calculate highest posterior density (HPD) of array for given credible_interval.
    The HPD is the minimum width Bayesian credible interval (BCI). This implementation works only
    for unimodal distributions.

    :param ary : posterior samples
    :param confidence_level : confidence level

    :return: intervals minima, intervals maxima
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
    """
    Computes properties of distribution based on its samples

    :param samples: samples of shape (n_samples, n_features)
    :param credible_intervals_levels: Confidence in (0, 1)
    of credible intervals to be computed
    :return: properties of distribution
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
    """
    Saves multi-clusters DE in an xlsx sheet

    :param filepath: xslx save path
    :param de_results: list of pandas Dataframes for each cluster
    :param cluster_names: list of cluster names
    :return:
    """
    writer = pd.ExcelWriter(filepath, engine="xlsxwriter")
    for i, x in enumerate(cluster_names):
        de_results[i].to_excel(writer, sheet_name=str(x))
    writer.close()
