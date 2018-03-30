import numpy as np 
from sklearn.neighbors import NearestNeighbors
import scipy.sparse
from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans
from sklearn.metrics import normalized_mutual_info_score as NMI
from sklearn.metrics import adjusted_rand_score as ARI
from sklearn.linear_model import LinearRegression
from sklearn.metrics import f1_score
from sklearn.metrics import roc_auc_score


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
    return score / 50.

def cluster_scores(latent_space, K, labels_true):
    labels_pred = KMeans(K, n_jobs=8, n_init=200).fit_predict(latent_space)
    return [silhouette_score(latent_space, labels_true), NMI(labels_true, labels_pred), ARI(labels_true, labels_pred)]



def dropout(X, rate=0.1):
    """
    X: original testing set
    ========
    returns:
    X_zero: copy of X with zeros
    i, j, ix: indices of where dropout is applied
    """
    X_zero = np.copy(X)
    # select non-zero subset
    i,j = np.nonzero(X_zero)
    
    # choice number 1 : select 10 percent of the non zero values (so that distributions overlap enough)
    ix = np.random.choice(range(len(i)), int(np.floor(0.1 * len(i))), replace=False)
    X_zero[i[ix], j[ix]] *= np.random.binomial(1, rate)
       
    # choice number 2, focus on a few but corrupt binomially
    #ix = np.random.choice(range(len(i)), int(slice_prop * np.floor(len(i))), replace=False)
    #X_zero[i[ix], j[ix]] = np.random.binomial(X_zero[i[ix], j[ix]].astype(np.int), rate)
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
    return np.median(np.abs(x - y))


# Evaluate AUC score DE
def auc_score_threshold(original_p_value, estimated_score, rank, p_value=True):
    # put ones on the smallest p_values
    true_labels = np.zeros_like(original_p_value)
    true_labels[original_p_value.argsort()[:rank]] = 1
    
    # make sure we are looking at the good ranking and not its inverse
    if p_value:
        estimated_score = -np.log(estimated_score)
    indices = np.isfinite(estimated_score)
    return roc_auc_score(true_labels[indices], estimated_score[indices])   
    
