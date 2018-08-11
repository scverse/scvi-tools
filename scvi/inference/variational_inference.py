import copy

import matplotlib.pyplot as plt
import numpy as np
import scipy
import torch
from scipy.stats import kde
from sklearn import neighbors
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.metrics import adjusted_rand_score as ARI
from sklearn.metrics import normalized_mutual_info_score as NMI
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture as GMM
from sklearn.neighbors import NearestNeighbors
from sklearn.utils.linear_assignment_ import linear_assignment
from torch.nn import functional as F

from scvi.dataset import CortexDataset
from scvi.dataset import GeneExpressionDataset
from scvi.dataset.data_loaders import DataLoaderWrapper
from scvi.dataset.data_loaders import TrainTestDataLoaders, TrainTestDataLoadersFish
from scvi.models.log_likelihood import compute_log_likelihood, compute_marginal_log_likelihood
from . import Inference

plt.switch_backend('agg')

from scipy.stats import itemfreq, entropy


def entropy_batch_mixing(latent_space, batches, n_neighbors=50, n_pools=50, n_samples_per_pool=100):
    batches = batches.ravel()
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(latent_space)
    indices = nbrs.kneighbors(latent_space, return_distance=False)[:, 1:]
    batch_indices = np.vectorize(lambda i: batches[i])(indices)
    entropies = np.apply_along_axis(entropy_from_indices, axis=1, arr=batch_indices)

    # average n_pools entropy results where each result is an average of n_samples_per_pool random samples.
    if n_pools == 1:
        score = np.mean(entropies)
    else:
        score = np.mean([
            np.mean(entropies[np.random.choice(len(entropies), size=n_samples_per_pool)])
            for _ in range(n_pools)
        ])

    return score


def entropy_from_indices(indices):
    return entropy(np.array(itemfreq(indices)[:, 1].astype(np.int32)))


def knn_purity(latent, label, n_neighbors=30):
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(latent)
    indices = nbrs.kneighbors(latent, return_distance=False)[:, 1:]
    neighbors_labels = np.vectorize(lambda i: label[i])(indices)

    # pre cell purity scores
    scores = ((neighbors_labels - label.reshape(-1, 1)) == 0).mean(axis=1)
    res = [np.mean(scores[label == i]) for i in np.unique(label)]  # per cell-type purity

    return np.mean(res)

def nn_overlap(X1, X2, k=100):
    nne = NearestNeighbors(n_neighbors=k + 1, n_jobs=8)
    assert len(X1) == len(X2)
    n_samples = len(X1)
    nne.fit(X1)
    kmatrix_1 = nne.kneighbors_graph(X1) - scipy.sparse.identity(n_samples)
    nne.fit(X2)
    kmatrix_2 = nne.kneighbors_graph(X2) - scipy.sparse.identity(n_samples)

    # 1 - spearman correlation from knn graphs
    spearman_correlation = scipy.stats.spearmanr(kmatrix_1.A.flatten(), kmatrix_2.A.flatten())[0]
    # 2 - fold enrichment
    set_1 = set(np.where(kmatrix_1.A.flatten() == 1)[0])
    set_2 = set(np.where(kmatrix_2.A.flatten() == 1)[0])
    fold_enrichment = len(set_1.intersection(set_2)) * n_samples ** 2 / (float(len(set_1)) * len(set_2))
    return spearman_correlation, fold_enrichment

def de_cortex(px_scale, all_labels, gene_names, M_permutation=100000, permutation=False):
    """
    Output average over statistics in a symmetric way (a against b)
    forget the sets if permutation is True
    :param M_permutation: 10000 - default value in Romain's code
    :param permutation:
    :return: A 1-d vector of statistics of size n_genes
    """
    # Compute sample rate for the whole dataset ?
    cell_types = np.array(['astrocytes_ependymal', 'endothelial-mural', 'interneurons', 'microglia',
                           'oligodendrocytes', 'pyramidal CA1', 'pyramidal SS'], dtype=np.str)
    # oligodendrocytes (#4) VS pyramidal CA1 (#5)
    couple_celltypes = (4, 5)  # the couple types on which to study DE

    print("\nDifferential Expression A/B for cell types\nA: %s\nB: %s\n" %
          tuple((cell_types[couple_celltypes[i]] for i in [0, 1])))

    # Here instead of A, B = 200, 400: we do on whole dataset then select cells
    sample_rate_a = (px_scale[all_labels.view(-1) == couple_celltypes[0]].view(-1, px_scale.size(1))
                     .cpu().detach().numpy())
    sample_rate_b = (px_scale[all_labels.view(-1) == couple_celltypes[1]].view(-1, px_scale.size(1))
                     .cpu().detach().numpy())

    # agregate dataset
    samples = np.vstack((sample_rate_a, sample_rate_b))

    # prepare the pairs for sampling
    list_1 = list(np.arange(sample_rate_a.shape[0]))
    list_2 = list(sample_rate_a.shape[0] + np.arange(sample_rate_b.shape[0]))
    if not permutation:
        # case1: no permutation, sample from A and then from B
        u, v = np.random.choice(list_1, size=M_permutation), np.random.choice(list_2, size=M_permutation)
    else:
        # case2: permutation, sample from A+B twice
        u, v = (np.random.choice(list_1 + list_2, size=M_permutation),
                np.random.choice(list_1 + list_2, size=M_permutation))

    # then constitutes the pairs
    first_set = samples[u]
    second_set = samples[v]

    res = np.mean(first_set >= second_set, 0)
    res = np.log(res + 1e-8) - np.log(1 - res + 1e-8)

    genes_of_interest = np.char.upper(["Thy1", "Mbp"])
    result = [(gene_name, res[np.where(gene_names == gene_name)[0]][0]) for gene_name in genes_of_interest]
    print('\n'.join([gene_name + " : " + str(r) for (gene_name, r) in result]))
    return result[1][1]  # if we had to give a metric to optimize


def plot_imputation(original, imputed, title="Imputation"):
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
    xi, yi = np.mgrid[0:ymax:nbins * 1j, 0:ymax:nbins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    plt.title(title, fontsize=12)
    plt.ylabel("Imputed counts")
    plt.xlabel('Original counts')

    plt.pcolormesh(yi, xi, zi.reshape(xi.shape), cmap="Reds")

    a, _, _, _ = np.linalg.lstsq(y[:, np.newaxis], x, rcond=-1)
    linspace = np.linspace(0, ymax)
    plt.plot(linspace, a * linspace, color='black')

    plt.plot(linspace, linspace, color='black', linestyle=":")
    plt.show()
    plt.savefig(title + '.png')


def proximity_imputation(real_latent1, normed_gene_exp_1, real_latent2, k=4):
    knn = neighbors.KNeighborsRegressor(k, weights='distance')
    y = knn.fit(real_latent1, normed_gene_exp_1).predict(real_latent2)
    return y


class VariationalInference(Inference):
    r"""The VariationalInference class for the unsupervised training of an autoencoder.

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SVAEC``
        :gene_dataset: A gene_dataset instance like ``CortexDataset()``
        :train_size: The train size, either a float between 0 and 1 or and integer for the number of training samples
         to use Default: ``0.8``.
        :\*\*kwargs: Other keywords arguments from the general Inference class.

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
        ... n_labels=gene_dataset.n_labels)

        >>> infer = VariationalInference(gene_dataset, vae, train_size=0.5)
        >>> infer.train(n_epochs=20, lr=1e-3)
    """
    default_metrics_to_monitor = ['ll']

    def __init__(self, model, gene_dataset, train_size=0.8, **kwargs):
        super(VariationalInference, self).__init__(model, gene_dataset, **kwargs)
        self.kl = None
        self.data_loaders = TrainTestDataLoaders(self.gene_dataset, train_size=train_size, use_cuda=self.use_cuda)

    def loss(self, tensors):
        sample_batch, local_l_mean, local_l_var, batch_index, _ = tensors
        reconst_loss, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var, batch_index)
        loss = torch.mean(reconst_loss + self.kl_weight * kl_divergence)
        return loss

    def on_epoch_begin(self):
        self.kl_weight = self.kl if self.kl is not None else min(1, self.epoch / self.n_epochs)

    def ll(self, name, verbose=False):
        ll = compute_log_likelihood(self.model, self.data_loaders[name])
        if verbose:
            print("LL for %s is : %.4f" % (name, ll))
        return ll

    ll.mode = 'min'

    def marginal_ll(self, name, verbose=False, n_mc_samples=1000):
        ll = compute_marginal_log_likelihood(self.model, self.data_loaders[name], n_mc_samples)
        if verbose:
            print("True LL for %s is : %.4f" % (name, ll))
        return ll

    def imputation_task(self, name, rate=0.1, n_samples=1, n_epochs=1, corruption="uniform"):
        corrupted_data = copy.deepcopy(self.gene_dataset.X)

        if corruption == "uniform":  # multiply the entry n with a Ber(0.9) random variable.
            i, j = np.nonzero(corrupted_data)
            ix = np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False)
            i, j = i[ix], j[ix]
            corrupted_data[i, j] *= np.random.binomial(n=np.ones(len(ix), dtype=np.int64), p=0.9)
        elif corruption == "binomial":  # multiply the entry n with a Bin(n, 0.9) random variable.
            i, j = (k.ravel() for k in np.indices(corrupted_data.shape))
            ix = np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False)
            i, j = i[ix], j[ix]
            corrupted_data[i, j] = np.random.binomial(n=corrupted_data[i, j].astype(np.int64), p=0.2)

        self.gene_dataset = gene_dataset = GeneExpressionDataset(
            *GeneExpressionDataset.get_attributes_from_matrix(
                corrupted_data,
                batch_indices=self.gene_dataset.batch_indices,
                labels=self.gene_dataset.labels
            )
        )

        original_data_loaders_loop = self.data_loaders.loop
        self.data_loaders.loop = ['corrupted_%s' % s for s in self.data_loaders.loop]
        original_keys = list(self.data_loaders.dict.keys())
        for key in original_keys:
            kwargs = copy.copy(self.data_loaders.kwargs)
            kwargs['collate_fn'] = gene_dataset.collate_fn
            kwargs['sampler'] = copy.copy(self.data_loaders[key].sampler)
            self.data_loaders['corrupted_%s' % key] = DataLoaderWrapper(gene_dataset, use_cuda=self.use_cuda, **kwargs)

        self.train(n_epochs=n_epochs)
        self.data_loaders.loop = original_data_loaders_loop

        original_list = []
        imputed_list = []
        batch_size = self.data_loaders.kwargs["batch_size"] // n_samples
        for tensors, corrupted_tensors in \
            zip(self.data_loaders[name].sequential(batch_size=batch_size),
                self.data_loaders['corrupted_%s' % name].sequential(batch_size=batch_size)):
            batch = tensors[0]
            actual_batch_size = batch.size(0)
            dropout_batch, _, _, batch_index, labels = corrupted_tensors
            px_rate = self.model.get_sample_rate(dropout_batch, batch_index=batch_index, y=labels, n_samples=n_samples)

            indices_dropout = torch.nonzero(batch - dropout_batch)
            i = indices_dropout[:, 0]
            j = indices_dropout[:, 1]

            batch = batch.unsqueeze(0).expand((n_samples, batch.size(0), batch.size(1)))
            original = np.array(batch[:, i, j].view(-1).cpu())
            imputed = np.array(px_rate[:, i, j].view(-1).cpu())

            cells_index = np.tile(np.array(i.cpu()), n_samples)

            original_list += [original[cells_index == i] for i in range(actual_batch_size)]
            imputed_list += [imputed[cells_index == i] for i in range(actual_batch_size)]

        return original_list, imputed_list

    def imputation(self, name, verbose=False, rate=0.1, n_samples=1, n_epochs=1, corruption="uniform",
                   title="Imputation"):
        original_list, imputed_list = self.imputation_task(name, rate=rate, n_epochs=n_epochs,
                                                           n_samples=n_samples, corruption=corruption)
        # Median of medians for all distances
        median_imputation_score = np.median(np.abs(np.concatenate(original_list) - np.concatenate(imputed_list)))

        # Mean of medians for each cell
        imputation_cells = []
        for original, imputed in zip(original_list, imputed_list):
            has_imputation = len(original) and len(imputed)
            imputation_cells += [np.median(np.abs(original - imputed)) if has_imputation else 0]
        mean_imputation_score = np.mean(imputation_cells)

        if verbose:
            print("Imputation Scores [corruption:%s - rate:%.2f] on  %s after %i:"
                  "\nMedian of Median: %.4f\nMean of Median for each cell: %.4f" %
                  (corruption, rate, name, n_epochs, median_imputation_score, mean_imputation_score))

        plot_imputation(np.concatenate(original_list), np.concatenate(imputed_list), title=title)
        return original_list, imputed_list

    def clustering_scores(self, name, verbose=True, prediction_algorithm='knn'):
        if self.gene_dataset.n_labels > 1:
            latent, _, labels = self.get_latent(name)
            if prediction_algorithm == 'knn':
                labels_pred = KMeans(self.gene_dataset.n_labels, n_init=200).fit_predict(latent)  # n_jobs>1 ?
            elif prediction_algorithm == 'gmm':
                gmm = GMM(self.gene_dataset.n_labels)
                gmm.fit(latent)
                labels_pred = gmm.predict(latent)

            asw_score = silhouette_score(latent, labels)
            nmi_score = NMI(labels, labels_pred)
            ari_score = ARI(labels, labels_pred)
            uca_score = unsupervised_clustering_accuracy(labels, labels_pred)[0]
            if verbose:
                print("Clustering Scores for %s:\nSilhouette: %.4f\nNMI: %.4f\nARI: %.4f\nUCA: %.4f" %
                      (name, asw_score, nmi_score, ari_score, uca_score))
            return asw_score, nmi_score, ari_score, uca_score

    def nn_overlap_score(self, name='sequential', verbose=True, **kwargs):
        if hasattr(self.gene_dataset, 'adt_expression_clr'):
            assert name == 'sequential'  # only works for the sequential data_loader (mapping indices)
            latent, _, _ = self.get_latent(name)
            protein_data = self.gene_dataset.adt_expression_clr
            spearman_correlation, fold_enrichment = nn_overlap(latent, protein_data, **kwargs)
            if verbose:
                print("Overlap Scores for %s:\nSpearman Correlation: %.4f\nFold Enrichment: %.4f" %
                      (name, spearman_correlation, fold_enrichment))
            return spearman_correlation, fold_enrichment

    def differential_expression_stats(self, name, *args, **kwargs):
        return self.de_stats(name, *args, **kwargs)

    def differential_expression(self, name, *args, verbose=False, **kwargs):
        px_scale, all_labels = self.differential_expression_stats(name, *args, **kwargs)

        if type(self.gene_dataset) == CortexDataset:
            if 'use_cuda' in kwargs:
                kwargs.pop('use_cuda')
            de_score = de_cortex(px_scale, all_labels, self.gene_dataset.gene_names, **kwargs)
            if verbose:
                print("DE score for cortex on %s is : %.4f" % (name, de_score))
        return de_score

    differential_expression.mode = 'max'

    def entropy_batch_mixing(self, name, verbose=False, **kwargs):
        if self.gene_dataset.n_batches == 2:
            latent, batch_indices, labels = self.get_latent(name)
            be_score = entropy_batch_mixing(latent, batch_indices, **kwargs)
            if verbose:
                print("Entropy batch mixing :", be_score)
            return be_score

    entropy_batch_mixing.mode = 'max'

    def knn_purity(self, name, verbose=False):
        latent, _, labels = get_latent(self.model, self.data_loaders[name])
        score = knn_purity(latent, labels)
        if verbose:
            print("KNN purity score :", score)
        return score

    knn_purity.mode = 'max'

    def adapt_encoder(self, n_path=10, n_epochs=50, frequency=5):
        vae = self.model
        params = list(vae.z_encoder.parameters()) + list(vae.l_encoder.parameters())
        z_encoder_state = copy.deepcopy(vae.z_encoder.state_dict())
        l_encoder_state = copy.deepcopy(vae.l_encoder.state_dict())
        self.data_loaders.loop = ['test']
        self.data_loaders.to_monitor = ['test']
        self.frequency = frequency

        # Training the model
        for i in range(n_path):
            # Re-initialize to create new path
            vae.z_encoder.load_state_dict(z_encoder_state)
            vae.l_encoder.load_state_dict(l_encoder_state)
            self.train(n_epochs, params=params)
        return min(self.history["ll_test"])

    def get_latent(self, name):
        latent = []
        batch_indices = []
        labels = []
        for tensors in self.data_loaders[name]:
            sample_batch, local_l_mean, local_l_var, batch_index, label = tensors
            latent += [self.model.sample_from_posterior_z(sample_batch, y=label)]
            batch_indices += [batch_index]
            labels += [label]
        return np.array(torch.cat(latent)), np.array(torch.cat(batch_indices)), np.array(torch.cat(labels)).ravel()

    def de_stats(self, name, M_sampling=100):
        """
        Output average over statistics in a symmetric way (a against b)
        forget the sets if permutation is True
        :param vae: The generative vae and encoder network
        :param data_loader: a data loader for a particular dataset
        :param M_sampling: number of samples
        :return: A 1-d vector of statistics of size n_genes
        """
        px_scales = []
        all_labels = []
        for tensors in self.data_loaders[name]:
            sample_batch, _, _, batch_index, labels = tensors
            sample_batch = sample_batch.repeat(1, M_sampling).view(-1, sample_batch.size(1))
            batch_index = batch_index.repeat(1, M_sampling).view(-1, 1)
            labels = labels.repeat(1, M_sampling).view(-1, 1)
            px_scales += [
                (self.model.get_sample_scale(sample_batch, batch_index=batch_index, y=labels).squeeze()).cpu()]
            all_labels += [labels.cpu()]

        px_scale = torch.cat(px_scales)
        all_labels = torch.cat(all_labels)

        return px_scale, all_labels

    def show_t_sne(self, name, n_samples=1000, color_by='', save_name='', latent=None, batch_indices=None,
                   labels=None, n_batch=None):
        # If no latent representation is given
        if latent is None:
            latent, batch_indices, labels = self.get_latent(name)
            latent, idx_t_sne = self.apply_t_sne(latent, n_samples)
            batch_indices = batch_indices[idx_t_sne].ravel()
            labels = labels[idx_t_sne].ravel()
        if not color_by:
            plt.figure(figsize=(10, 10))
            plt.scatter(latent[:, 0], latent[:, 1])
        if color_by == 'scalar':
            plt.figure(figsize=(10, 10))
            plt.scatter(latent[:, 0], latent[:, 1], c=labels.ravel())
        else:
            if n_batch is None:
                n_batch = self.gene_dataset.n_batches
            if color_by == 'batches' or color_by == 'labels':
                indices = batch_indices if color_by == 'batches' else labels
                n = n_batch if color_by == 'batches' else self.gene_dataset.n_labels
                if hasattr(self.gene_dataset, 'cell_types') and color_by == 'labels':
                    plt_labels = self.gene_dataset.cell_types
                else:
                    plt_labels = [str(i) for i in range(len(np.unique(indices)))]
                plt.figure(figsize=(10, 10))
                for i, label in zip(range(n), plt_labels):
                    plt.scatter(latent[indices == i, 0], latent[indices == i, 1], label=label)
                plt.legend()
            elif color_by == 'batches and labels':
                fig, axes = plt.subplots(1, 2, figsize=(14, 7))
                for i in range(n_batch):
                    axes[0].scatter(latent[batch_indices == i, 0], latent[batch_indices == i, 1], label=str(i))
                axes[0].set_title("batch coloring")
                axes[0].axis("off")
                axes[0].legend()

                indices = labels
                if hasattr(self.gene_dataset, 'cell_types'):
                    plt_labels = self.gene_dataset.cell_types
                else:
                    plt_labels = [str(i) for i in range(len(np.unique(indices)))]
                for i, cell_type in zip(range(self.gene_dataset.n_labels), plt_labels):
                    axes[1].scatter(latent[indices == i, 0], latent[indices == i, 1], label=cell_type)
                axes[1].set_title("label coloring")
                axes[1].axis("off")
                axes[1].legend()
        plt.axis("off")
        plt.tight_layout()
        if save_name:
            plt.savefig(save_name)

    @staticmethod
    def apply_t_sne(latent, n_samples=1000):
        idx_t_sne = np.random.permutation(len(latent))[:n_samples] if n_samples else np.arange(len(latent))
        if latent.shape[1] != 2:
            latent = TSNE().fit_transform(latent[idx_t_sne])
        return latent, idx_t_sne


class VariationalInferenceFish(VariationalInference):
    r"""The VariationalInference class for the unsupervised training of an autoencoder.

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SVAEC``
        :gene_dataset: A gene_dataset instance like ``CortexDataset()``
        :train_size: The train size, either a float between 0 and 1 or and integer for the number of training samples
         to use Default: ``0.8``.
        :\*\*kwargs: Other keywords arguments from the general Inference class.

    Examples:
        >>> gene_dataset_seq = CortexDataset()
        >>> gene_dataset_fish = SmfishDataset()
        >>> vae = VAE(gene_dataset_seq.nb_genes, gene_dataset_fish.nb_genes,
        ... n_labels=gene_dataset.n_labels, use_cuda=True)

        >>> infer = VariationalInference(gene_dataset_seq, gene_dataset_fish, vae, train_size=0.5)
        >>> infer.train(n_epochs=20, lr=1e-3)
    """
    default_metrics_to_monitor = ['ll']

    def __init__(self, model, gene_dataset_seq, gene_dataset_fish, train_size=0.8, use_cuda=True, cl_ratio=0,
                 n_epochs_even=1, n_epochs_kl=2000, n_epochs_cl=1, **kwargs):
        super(VariationalInferenceFish, self).__init__(model, gene_dataset_seq, use_cuda=use_cuda, **kwargs)
        self.kl = None
        self.cl_ratio = cl_ratio
        self.n_epochs_cl = n_epochs_cl
        self.n_epochs_even = n_epochs_even
        self.n_epochs_kl = n_epochs_kl
        self.weighting = 0
        self.kl_weight = 0
        self.classification_ponderation = 0
        self.data_loaders = TrainTestDataLoadersFish(gene_dataset_seq, gene_dataset_fish,
                                                     train_size=train_size, use_cuda=use_cuda)

    def loss(self, tensors_seq, tensors_fish):
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors_seq
        reconst_loss, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var, batch_index, mode="scRNA",
                                                 weighting=self.weighting)
        # If we want to add a classification loss
        if self.cl_ratio != 0:
            reconst_loss += self.cl_ratio * F.cross_entropy(self.model.classify(sample_batch, mode="scRNA"),
                                                            labels.view(-1))
        loss = torch.mean(reconst_loss + self.kl_weight * kl_divergence)
        sample_batch_fish, local_l_mean, local_l_var, batch_index_fish, _, _, _ = tensors_fish
        reconst_loss_fish, kl_divergence_fish = self.model(sample_batch_fish, local_l_mean, local_l_var,
                                                           batch_index_fish, mode="smFISH")
        loss_fish = torch.mean(reconst_loss_fish + self.kl_weight * kl_divergence_fish)
        loss = loss * sample_batch.size(0) + loss_fish * sample_batch_fish.size(0)
        loss /= (sample_batch.size(0) + sample_batch_fish.size(0))
        return loss + loss_fish

    def on_epoch_begin(self):
        self.weighting = min(1, self.epoch / self.n_epochs_even)
        self.kl_weight = self.kl if self.kl is not None else min(1, self.epoch / self.n_epochs_kl)
        self.classification_ponderation = min(1, self.epoch / self.n_epochs_cl)

    def ll(self, name, verbose=False):
        data_loader = self.data_loaders[name]
        log_lkl = 0
        for i_batch, tensors in enumerate(data_loader):
            if name == "train_fish" or name == "test_fish":
                sample_batch, local_l_mean, local_l_var, batch_index, labels, _, _ = tensors
                reconst_loss, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var,
                                                         batch_index=batch_index,
                                                         y=labels, mode="smFISH")
            else:
                sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors
                reconst_loss, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var,
                                                         batch_index=batch_index, y=labels)
            log_lkl += torch.sum(reconst_loss).item()
        n_samples = (len(data_loader.dataset)
                     if not (hasattr(data_loader, 'sampler') and hasattr(data_loader.sampler, 'indices')) else
                     len(data_loader.sampler.indices))
        ll = log_lkl / n_samples
        if verbose:
            print("LL for %s is : %.4f" % (name, ll))
        return ll

    def show_spatial_expression(self, x_coord, y_coord, labels, color_by='scalar', title='spatial_expression.svg'):
        x_coord = x_coord.reshape(-1, 1)
        y_coord = y_coord.reshape(-1, 1)
        latent = np.concatenate((x_coord, y_coord), axis=1)
        self.show_t_sne(name=None, n_samples=1000, color_by=color_by, save_name=title,
                        latent=latent, batch_indices=None,
                        labels=labels)


def unsupervised_clustering_accuracy(y, y_pred):
    """
    Unsupervised Clustering Accuracy
    """
    assert len(y_pred) == len(y)
    n_clusters = len(np.unique(y))
    reward_matrix = np.zeros((n_clusters, n_clusters), dtype=np.int64)
    for y_pred_, y_ in zip(y_pred, y):
        reward_matrix[y_pred_, y_] += 1
    cost_matrix = reward_matrix.max() - reward_matrix
    ind = linear_assignment(cost_matrix)
    return sum([reward_matrix[i, j] for i, j in ind]) * 1.0 / y_pred.size, ind
