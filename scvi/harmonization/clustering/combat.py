import numpy as np
import rpy2.robjects as ro

import warnings
from rpy2.rinterface import RRuntimeWarning
import rpy2.robjects.numpy2ri as numpy2ri
from scipy.io import mmwrite
from sklearn.decomposition import PCA

class COMBAT():
    def __init__(self):
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        numpy2ri.activate()
        ro.r["library"]("gmodels")
        ro.r["library"]("sva")
        ro.r["library"]("Matrix")
        ro.r["library"]("RcppCNPy")
        # ro.r["library"]("reticulate")
        # save_npz('temp.npz', csr_matrix([[1, 2, 0], [0, 0, 3], [4, 0, 5]]))
        # ro.r('sparse <- import("scipy.sparse")')
        # ro.r('X <- sparse$load_npz("temp.npz")')

    def csr2r(self, matrix):
        # because rpy2 don't have sparse encoding try printing it to mtx and reading it in R
        # the object is named X
        mmwrite('temp.mtx',matrix)
        ro.r('X <- readMM("temp.mtx")')
        # save_npz('temp.npz', matrix)
        # ro.r('sparse <- import("scipy.sparse")')
        # ro.r('X <- sparse$load_npz("temp.npz")')

    def combat_correct(self, dataset):
        batch_indices = np.concatenate(dataset.batch_indices)
        ro.r.assign("batch", ro.IntVector(batch_indices))
        # X = np.asarray(dataset.X.T.todense())
        # nb_genes, n_samples = X.shape
        # X = ro.r.matrix(X, nrow=nb_genes, ncol=n_samples)
        # ro.r.assign('X', X)
        self.csr2r(dataset.X.T)
        corrected = ro.r('ComBat(log(as.matrix(X)+1),batch)')
        return corrected

    def combat_pca(self, dataset):
        corrected = self.combat_correct(dataset)
        pca = PCA(n_components=10)
        pca.fit(corrected)
        pc = pca.components_
        return pc

