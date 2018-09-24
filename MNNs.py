import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import warnings
from rpy2.rinterface import RRuntimeWarning


class MNN(BaseEstimator, TransformerMixin):

    def __init__(self, n_components=10, learn_V=True):
        self.n_components = n_components
        self.learn_V = learn_V
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        rpy2.robjects.numpy2ri.activate()
        ro.r["library"]("scran")
        ro.r["library"]("BiocParallel")
        ro.r["library"]("scRNAseq")
        ro.r["library"]("matrixStats")
        ro.r["library"]("magrittr")
        ro.r["library"]("ggplot2")
        ro.r["library"]("biomaRt")
        ro.r["library"]("tibble")
        ro.r["library"]("SIMLR")
        ro.r("BiocParallel::register(BiocParallel::MulticoreParam(4))")

    def fit_transform(self, X, batch, list_b):
        index_0 = np.where(batch == list_b[0])[0]
        index_1 = np.where(batch == list_b[1])[0]

        self.A_ = np.log(1 + X[index_0].T)
        self.B_ = np.log(1 + X[index_1].T)

        nr, nc = self.A_.shape
        Ar = ro.r.matrix(self.A_, nrow=nr, ncol=nc)
        ro.r.assign("matrix_A", Ar)

        nr, nc = self.B_.shape
        Br = ro.r.matrix(self.B_, nrow=nr, ncol=nc)
        ro.r.assign("matrix_B", Br)

        ro.r("out <- mnnCorrect(matrix_A, matrix_B, BPPARAM=MulticoreParam(4))")

        corr_A = np.array(ro.r("out$corrected")[0]).T
        corr_B = np.array(ro.r("out$corrected")[1]).T

        arr = np.zeros_like(X, dtype=np.float)
        arr[index_0] = corr_A
        arr[index_1] = corr_B

        return arr
