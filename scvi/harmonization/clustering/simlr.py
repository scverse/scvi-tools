import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_array, check_is_fitted
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import warnings
from rpy2.rinterface import RRuntimeWarning

class SIMLR(BaseEstimator, TransformerMixin):
    # An algorithm for clustering (to compare to KNN in the latent space)

    def __init__(self, n_clusters=10):
        self.n_clusters = n_clusters
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        rpy2.robjects.numpy2ri.activate()
        ro.r["library"]("SIMLR")
        ro.r["library"]("BiocParallel")
        ro.r["library"]("matrixStats")
        ro.r["library"]("magrittr")
        ro.r["library"]("ggplot2")
        ro.r["library"]("biomaRt")
        ro.r["library"]("tibble")
        ro.r("BiocParallel::register(BiocParallel::MulticoreParam())")
        ro.r.assign("K", n_clusters)


    def estimate_clusters_numbers(self, X):
        nr,nc = X.shape
        X_trainr = ro.r.matrix(X, nrow=nr, ncol=nc)
        ro.r.assign("matrix_val", X_trainr)
        ro.r("NUMC = 2:20")
        ro.r("out <- SIMLR_Estimate_Number_of_Clusters(t(log(1+matrix_val)), NUMC = NUMC, cores.ratio = 0)")
        self.estimated_clusters = ro.r("out$K1")

    def fit_transform(self, X):
        self.X_ = X
        nr,nc = X.shape
        X_trainr = ro.r.matrix(X, nrow=nr, ncol=nc)
        ro.r.assign("matrix_val", X_trainr)

        ro.r("out <- SIMLR_Large_Scale(t(log(1+matrix_val)), c = K, k=30,kk=200)")
        self.S = ro.r("out$S")
        self.S0 = ro.r("out$S0")
        self.F = ro.r("out$F")

        self.clusters = ro.r("out$y$cluster")
        self.ydata = ro.r("out$ydata")
        # Return the classifier
        return self


    def output_estimation(self):
        """
        Returns parameters
        """
        return self.params
