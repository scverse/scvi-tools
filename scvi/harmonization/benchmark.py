import numpy as np
from scvi.dataset import GeneExpressionDataset
from scvi.metrics.clustering import get_latent, entropy_batch_mixing

from scipy import sparse
from copy import deepcopy

def sample_by_batch(batch_indices, nsamples):
    nbatches = len(np.unique(batch_indices))
    if isinstance(nsamples,int):
        nsamples = np.repeat(nsamples, nbatches)
    sample = []
    for i in np.unique(batch_indices):
        idx = np.arange(len(batch_indices))[batch_indices == i]
        s = np.random.permutation(idx)[:min(len(idx), nsamples[i])]
        sample.append(s)
    sample = np.concatenate(sample)
    return sample

def harmonization_stat(model, data_loader,keys, pop1, pop2):
    latent, batch_indices, labels = get_latent(model, data_loader)
    batch_indices = np.concatenate(batch_indices)
    sample = sample_by_batch(batch_indices, 2000)
    sample_2batch = sample[(batch_indices[sample] == pop1) + (batch_indices[sample] == pop2)]
    batch_entropy = entropy_batch_mixing(latent[sample_2batch, :], batch_indices[sample_2batch])
    print("Entropy batch mixing : %f.3" % batch_entropy)
    sample = sample_by_batch(labels, 200)
    res = knn_purity_avg(
        latent[sample,:], labels.astype('int')[sample],
        keys, acc=True
    )
    print("Average knn purity : %f.3" % np.mean([x[1] for x in res]))
    return(batch_entropy,res)


def ind_log_likelihood(vae, data_loader):
    # Iterate once over the data_loader and computes the total log_likelihood
    log_lkl = []
    for i_batch, tensors in enumerate(data_loader):
        if vae.use_cuda:
            tensors = to_cuda(tensors)
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors
        sample_batch = sample_batch.type(torch.float32)
        reconst_loss, kl_divergence = vae(sample_batch, local_l_mean, local_l_var, batch_index=batch_index,
                                          y=labels)
        log_lkl.append(torch.sum(reconst_loss).item())
    return log_lkl


def assign_label(cellid, geneid, labels_map, count, cell_type, seurat):
    labels = seurat[1:, 4]
    labels = np.int64(np.asarray(labels))
    labels_new = deepcopy(labels)
    for i, j in enumerate(labels_map):
        labels_new[labels == i] = j
    temp = dict(zip(cellid, count))
    new_count = []
    for x in seurat[1:, 5]:
        new_count.append(temp[x])
    new_count = sparse.vstack(new_count)
    dataset = GeneExpressionDataset(*GeneExpressionDataset.get_attributes_from_matrix(new_count, labels=labels_new),
                                    gene_names=geneid, cell_types=cell_type)
    return dataset

def knn_purity(latent, label, batch_indices, pop1, pop2, keys, n_sample=1000,acc=False):
    sample = np.random.permutation(len(label))[range(n_sample)]
    latent = latent[sample]
    label = label[sample]
    batch_indices = batch_indices[sample]
    if str(type(latent)) == "<class 'scipy.sparse.csr.csr_matrix'>":
        latent = latent.todense()
    distance = np.zeros((n_sample, n_sample))
    neighbors_graph = np.zeros((n_sample, n_sample))
    for i in range(n_sample):
        for j in range(i, n_sample):
            distance[i, j] = distance[j, i] = np.sum(np.asarray(latent[i] - latent[j]) ** 2)
    for i, d in enumerate(distance):
        neighbors_graph[i, d.argsort()[:30]] = 1
    kmatrix = neighbors_graph - np.identity(latent.shape[0])
    score = []
    for i in range(n_sample):
        if batch_indices[i] == pop1:
            lab = label[i]
            n_lab = label[(kmatrix[i] == 1) * (batch_indices == pop2)]
            if len(n_lab) > 0:
                if acc is False:
                    score.append(np.sum([x == lab for x in n_lab]) / len(n_lab))
                else:
                    if np.sum([x == lab for x in n_lab]) / len(n_lab) > 0.5:
                        score.append(1)
                    else:
                        score.append(0)
            else:
                score.append(-1)
    score = np.asarray(score)
    label = label[batch_indices == pop1]
    label = label[score != (-1)]
    score = score[score != (-1)]
    keys = keys[np.unique(label[batch_indices is pop1])]
    res = [np.mean(np.asarray(score)[label == i]) for i in np.unique(label)]
    res = [[keys[i], res[i]] for i in range(len(res))]
    return res


def knn_purity_avg(latent, label, keys, n_sample=1000,acc=False):
    sample = np.random.permutation(len(label))[range(n_sample)]
    latent = latent[sample]
    label = label[sample]
    if str(type(latent)) == "<class 'scipy.sparse.csr.csr_matrix'>":
        latent = latent.todense()
    distance = np.zeros((n_sample, n_sample))
    neighbors_graph = np.zeros((n_sample, n_sample))
    for i in range(n_sample):
        for j in range(i, n_sample):
            distance[i, j] = distance[j, i] = np.sum(np.asarray(latent[i] - latent[j]) ** 2)
    for i, d in enumerate(distance):
        neighbors_graph[i, d.argsort()[:30]] = 1
    kmatrix = neighbors_graph - np.identity(latent.shape[0])
    score = []
    for i in range(n_sample):
        lab = label[i]
        n_lab = label[kmatrix[i] == 1]
        if acc is False:
            score.append(np.sum([x == lab for x in n_lab]) / len(n_lab))
        else:
            if np.sum([x == lab for x in n_lab]) / len(n_lab) > 0.5:
                score.append(1)
            else:
                score.append(0)
    score = np.asarray(score)
    res = [np.mean(np.asarray(score)[label == i]) for i in np.unique(label)]
    res = [[keys[i], res[i], np.sum(label == k)] for i,k in enumerate(np.unique(label))]
    return res


def sample_celltype(dataset1, cell, prop):
    genedataset = deepcopy(dataset1)
    celltype = [i for i, name in enumerate(genedataset.cell_types) if name == cell]
    labs = np.concatenate(genedataset.labels)
    idx1 = np.arange(len(labs))[labs != celltype]
    idx2 = np.arange(len(labs))[labs == celltype]
    idx2 = np.random.permutation(idx2)[:int(len(idx2) * prop)]
    idx = np.concatenate([idx1, idx2])
    genedataset.X = genedataset.X[idx, :]
    genedataset.labels = genedataset.labels[idx]
    genedataset.batch_indices = genedataset.batch_indices[idx]
    return genedataset

