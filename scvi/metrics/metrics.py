import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.stats import itemfreq, entropy
use_cuda = True

def entropy_batch_mixing(latent_space, batches, n_neighbors=50):
    batches = batches.ravel()
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(latent_space)
    indices = nbrs.kneighbors(latent_space, return_distance=False)[:, 1:]
    batch_indices = np.vectorize(lambda i: batches[i])(indices)
    return np.mean(np.apply_along_axis(entropy_from_indices, axis=1, arr=batch_indices))

def entropy_from_indices(indices):
    return entropy(np.array(itemfreq(indices)[:, 1].astype(np.int32)))


def knn_purity_avg(latent, label, keys, n_sample=1000,acc=False):
    sample = np.random.permutation(len(label))[range(n_sample)]
    latent = latent[sample]
    label = label[sample]
    if str(type(latent)) == "<class 'scipy.sparse.csr.csr_matrix'>":
        latent = latent.todense()
    distance = np.zeros((n_sample, n_sample))


    n_neighbors = 30
    # neighbors_graph = np.zeros((n_sample, n_sample))
    # for i in range(n_sample):
    #     for j in range(i, n_sample):
    #         distance[i, j] = distance[j, i] = np.sum(np.asarray(latent[i] - latent[j]) ** 2)
    # for i, d in enumerate(distance):
    #     neighbors_graph[i, d.argsort()[:30]] = 1
    # kmatrix = neighbors_graph - np.identity(latent.shape[0])
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(latent)
    indices = nbrs.kneighbors(latent, return_distance=False)[:, 1:]
    labels = np.vectorize(lambda i: label[i])(indices)

    score = []
    for i in range(n_sample):
        lab = label[i]
        n_lab = labels[i]
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

# sample = np.random.permutation(len(label))[range(n_sample)]
# latent = latent[sample]
# label = label[sample]



def metrics_be_knn(latent,batch_indices,tm_facs_meta,tm_droplet_meta):
    n_samples=2000
    np.random.seed(1)
    indices = []
    for i in np.unique(batch_indices):
        indices_i = np.where(batch_indices.ravel() == i)[0]
        indices += [indices_i[np.random.permutation(len(indices_i))][:n_samples]]

    sample = np.concatenate(indices, axis=0)

    print("Entropy batch mixing :", entropy_batch_mixing(latent[sample], batch_indices[sample]))


    labels = np.concatenate([tm_facs_meta[tm_facs_meta.tissue=='Marrow'].cell_ontology_class,tm_droplet_meta[tm_droplet_meta.tissue=='Marrow'].cell_ontology_class])
    labels = [str(x) for x in labels]
    keys, labels = np.unique(labels, return_inverse=True)
    labels = np.asarray(labels)

    np.random.seed(1)
    n_samples= 1000
    indices = []
    for i in np.unique(labels):
        indices_i = np.where(labels == i)[0]
        indices += [indices_i[np.random.permutation(len(indices_i))][:n_samples]]

    sample = np.concatenate(indices, axis=0)


    res = knn_purity_avg(
        latent[sample, :], labels[sample],
        keys=keys, acc=True
    )

    print("avg acc", np.mean([x[1] for x in res]))
