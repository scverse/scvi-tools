import numpy as np 
from sklearn.neighbors import NearestNeighbors
import scipy.sparse
from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans
from sklearn.metrics import normalized_mutual_info_score as NMI
from sklearn.metrics import adjusted_rand_score as ARI
from sklearn.linear_model import LinearRegression
from sklearn.metrics import f1_score


# CLUSTERING METRICS
def entropy_batch_mixing(latent_space, batches):
    def entropy(hist_data):
        n_batches = len(np.unique(hist_data))
        if n_batches > 2:
            raise ValueError("Should be only two clusters for this metric")
        frequency = np.mean(hist_data == 1)
        if frequency == 0 or frequency == 1:
            return 0
        return -frequency * np.log(frequency) - (1 - frequency) * np.log(1 - frequency)

    nne = NearestNeighbors(n_neighbors=51, n_jobs=8)
    nne.fit(latent_space)
    kmatrix = nne.kneighbors_graph(latent_space) - scipy.sparse.identity(latent_space.shape[0])

    score = 0
    for t in range(50):
        indices = np.random.choice(np.arange(latent_space.shape[0]), size=100)
        score += np.mean([entropy(batches[kmatrix[indices].nonzero()[1]\
                                 [kmatrix[indices].nonzero()[0] == i]]) for i in range(100)])
    return score / 50

def nne(latent_space, clusters):
    nne = NearestNeighbors(n_neighbors=2, n_jobs=8)
    nne.fit(latent_space)
    kmatrix = nne.kneighbors_graph(latent_space) - scipy.sparse.identity(latent_space.shape[0])
    l = np.where(kmatrix.A > 0)
    return np.mean(clusters[l[0]] == clusters[l[1]])

def distances_metrics(latents, labels_true):
    return {"nne":nne(latents, labels_true), "silhouette":silhouette_score(latents, labels_true)}

def cluster_latent(latent_space, K):
    return KMeans(K, n_jobs=8).fit_predict(latent_space)

def clustering_metrics(labels_pred, labels_true):
    return {"NMI":NMI(labels_true, labels_pred), "ARI":ARI(labels_true, labels_pred),\
            "F1":f1_score(labels_true, labels_pred, average='weighted')}

# REPEATED CLUSTERING
def best_clusters_score(projection, K, clusters, run=50):
    res = clustering_metrics(cluster_latent(projection, K), clusters)
    for i in range(run):
        current_res = clustering_metrics(cluster_latent(projection, K), clusters)
        for key, value in res.iteritems():
            if value < current_res[key]:
                res[key] = current_res[key]
    return res


#ZERO GENERATION

def dropout(X, decay=0, uniform=True):
    """
    X: original testing set
    decay: decay parameter estimated by ZIFA
    ========
    returns:
    X_zero: copy of X with zeros
    i, j, ix: indices of where dropout is applied
    """
    X_zero = np.copy(X)
    # compute dropout prob
    p = np.exp( - decay * np.log(1 + X_zero)**2)
    # select non-zero subset
    i,j = np.nonzero(X_zero)
    ix = np.random.choice(range(len(i)), int(np.floor(0.1 * len(i))), replace=False)
    if uniform == False:
        rate = 1-p[i[ix], j[ix]]
    else:
        rate = 0.1
    X_zero[i[ix], j[ix]] *= np.random.binomial(1, rate )
    return X_zero, i, j, ix


# IMPUTATION METRICS
def imputation_error(X_mean, X, X_zero, i, j, ix):
    """
    X_mean: imputed dataset
    X: original dataset
    X_zero: zeros dataset
    i, j, ix: indices of where dropout was applied
    ========
    returns:
    median L1 distance between datasets at indices given
    """
    all_index = i[ix], j[ix]
    x, y = X_mean[all_index], X[all_index]
    return np.median(np.abs(x[X_zero[all_index] == 0] - y[X_zero[all_index] == 0]))

def mask_nan_inf(data):
    return data[np.logical_and(np.isinf(data) == False, np.isnan(data)==False)]

def imputation_dropout(log_dropout, X_zero, i, j, ix):
    """
    log_dropout: imputed dropout log prob
    X_zero: zeroed dataset
    i, j, ix: indices of where dropout was applied
    ========
    returns:
    average log likelihood of being dropped out
    """

    X_slice = X_zero[i[ix], j[ix]]
    log_slice = log_dropout[i[ix], j[ix]]
    
    ent = np.zeros_like(log_slice)
    ent[X_slice == 0] = -log_slice[X_slice == 0]
    inverse_log = np.log1p(- np.exp(log_slice))
    ent[X_slice > 0] = -inverse_log[X_slice >0]
    median_entropy = np.median(mask_nan_inf(ent))
    mean_entropy = np.mean(mask_nan_inf(ent))

    return median_entropy, mean_entropy

# Evaluate correlation with RUV
def correlation_metrics(latent_space, RUV):
    # regress RUV
    RUV_regressor = LinearRegression(n_jobs=8)
    RUV_regressor.fit(latent_space, RUV)
    r_score = RUV_regressor.score(latent_space, RUV)
    return r_score


# Evaluate AUC score DE
from sklearn.metrics import roc_auc_score
def auc_score_threshold(original_p_value, estimated_p_value, rank, log=True):
    reject = -np.log(original_p_value)
    true_labels = np.zeros_like(original_p_value)
    true_labels[original_p_value.argsort()[:rank]] = 1
    indices = np.logical_not(np.logical_or(np.isnan(estimated_p_value), np.isinf(estimated_p_value)))
    if log:
        return roc_auc_score(true_labels[indices], -np.log(estimated_p_value[indices]))
    else:
        return roc_auc_score(true_labels[indices], estimated_p_value[indices])
    
    