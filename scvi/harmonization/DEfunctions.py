import numpy as np
import pandas as pd

use_cuda = True


def GetScale(trainer, n_batches, n_samples_per_batch, gene_dataset, selection, seed=0):
    scale = []
    np.random.seed(seed)
    for batch in range(n_batches):
        idx = np.random.choice(np.arange(len(gene_dataset))[selection], n_samples_per_batch)
        posterior = trainer.create_posterior(trainer.model, gene_dataset, indices=idx)
        scale.append(posterior.sequential().get_harmonized_scale(batch))
    scale = np.concatenate(scale)
    return scale


def GetScaleFromBatch(trainer, batchid, n_samples_per_batch, gene_dataset, selection, seed=0):
    scale = []
    np.random.seed(seed)
    for batch in batchid:
        idx = np.random.choice(np.arange(len(gene_dataset))[selection], n_samples_per_batch)
        posterior = trainer.create_posterior(trainer.model, gene_dataset, indices=idx)
        scale.append(posterior.sequential().get_harmonized_scale(batch))
    scale = np.concatenate(scale)
    return scale


def GetBayesFactor(scale1, scale2, m_permutation):
    res = np.repeat(0, scale1.shape[1])
    for _ in range(m_permutation):
        sample1 = scale1[np.random.choice(np.arange(scale1.shape[0]), 1), :].ravel()
        sample2 = scale2[np.random.choice(np.arange(scale2.shape[0]), 1), :].ravel()
        res = res + (sample1 > sample2)
    res = res / m_permutation
    res = np.log(res + 1e-8) - np.log(1 - res + 1e-8)
    return res


def ComputeDE(trainer, gene_dataset, idx1, idx2, batch1, batch2, nrep=1, n_samples_per_batch=1000,
              m_permutation=10000, seed=0):
    bayes = []
    np.random.seed(seed)
    seeds = np.random.choice(1000000, nrep)
    for i in np.arange(0, nrep):
        scale1 = GetScaleFromBatch(trainer, batch1, n_samples_per_batch, gene_dataset, idx1, seed=seeds[i])
        scale2 = GetScaleFromBatch(trainer, batch2, n_samples_per_batch, gene_dataset, idx2, seed=seeds[i])
        bf = GetBayesFactor(scale1, scale2, m_permutation)
        bayes.append(bf)
    bayes = np.asarray(bayes)
    return bayes


def RawCountsProperties(gene_dataset, idx1, idx2):
    mean1 = np.mean(gene_dataset.X[idx1, :], axis=0)
    mean2 = np.mean(gene_dataset.X[idx2, :], axis=0)
    nonz1 = np.mean(gene_dataset.X[idx1, :] != 0, axis=0)
    nonz2 = np.mean(gene_dataset.X[idx2, :] != 0, axis=0)
    return np.asarray(mean1).ravel(), np.asarray(mean2).ravel(), np.asarray(nonz1).ravel(), np.asarray(nonz2).ravel()


def One2AllBayes(trainer, gene_dataset, clusters, batchid, subset=None, mincells=10):
    if subset is None:
        subset = [0]
    de_res = []
    de_cluster = []
    cluster_id = np.unique(clusters[clusters >= 0])
    for i in cluster_id:
        if subset is not None:
            idx1 = (clusters == i)
            idx2 = (clusters != i)
        else:
            idx1 = (clusters == i) * subset
            idx2 = (clusters != i) * subset
        if np.sum(idx1) > mincells and np.sum(idx2) > mincells:
            de_cluster.append(i)
            bayes = ComputeDE(trainer, gene_dataset, idx1, idx2, np.unique(batchid), np.unique(batchid), nrep=2)
            mean1, mean2, nonz1, nonz2 = RawCountsProperties(gene_dataset, idx1, idx2)
            res = pd.DataFrame([bayes[0, :], bayes[1, :], mean1, mean2, nonz1, nonz2],
                               index=['bayes1', 'bayes2', 'mean1', 'mean2', 'nonz1', 'nonz2'],
                               columns=gene_dataset.gene_names).T
            res = res.loc[np.asarray([not x.startswith('RPS') for x in res.index])]
            res = res.loc[np.asarray([not x.startswith('RPL') for x in res.index])]
            res = res.sort_values(by=['bayes1'], ascending=False)
            de_res.append(res)
    return de_res, de_cluster


def WithinClusterDE(trainer, gene_dataset, clusters, group1, group2, batch1, batch2, mincells=10):
    nc = len(np.unique(clusters[clusters >= 0]))
    de_res = []
    vc = []
    for i in np.arange(0, nc):
        idx1 = (clusters == i) * group1
        idx2 = (clusters == i) * group2
        if np.sum(idx1) > mincells and np.sum(idx2) > mincells:
            vc.append(i)
            bayes = ComputeDE(trainer, gene_dataset, idx1, idx2, batch1, batch2, nrep=2)
            mean1, mean2, nonz1, nonz2 = RawCountsProperties(gene_dataset, idx1, idx2)
            res = pd.DataFrame([bayes[0, :], bayes[1, :], mean1, mean2, nonz1, nonz2],
                               index=['bayes1', 'bayes2', 'mean1', 'mean2', 'nonz1', 'nonz2'],
                               columns=gene_dataset.gene_names).T
            up = res.loc[(res['bayes1'] > 0)]
            up = up.sort_values(by=['bayes1'], ascending=False)
            down = res.loc[(res['bayes1'] < 0)]
            down = down.sort_values(by=['bayes1'])
            de_res.append([up, down])
    return de_res, vc
