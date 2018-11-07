import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import warnings
import numpy as np
from rpy2.rinterface import RRuntimeWarning


class IDR(object):

    def __init__(self):
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        rpy2.robjects.numpy2ri.activate()
        ro.r["library"]("idr")

    def fit(self, p_val_1, p_val_2, p_prior=0.1):
        p_1 = np.copy(p_val_1)
        p_2 = np.copy(p_val_2)
        indices = np.logical_and(np.isfinite(p_1), np.isfinite(p_2))

        nr, nc = p_val_1[indices, np.newaxis].shape
        p_val_1r = ro.r.matrix(p_val_1[indices, np.newaxis], nrow=nr, ncol=nc)
        ro.r.assign("p_val_1", p_val_1r)

        nr, nc = p_val_2[indices, np.newaxis].shape
        p_val_2r = ro.r.matrix(p_val_2[indices, np.newaxis], nrow=nr, ncol=nc)
        ro.r.assign("p_val_2", p_val_2r)

        ro.r("x <- cbind(p_val_1[, 1], p_val_2[, 1])")
        ro.r("mu = 1")
        ro.r("sigma = 0.5")
        ro.r("rho = 0.5")
        ro.r.assign("p", 0.1)
        ro.r("idr.out <- est.IDR(x, mu, sigma, rho, p, eps=0.001, max.ite=20)")
        return ro.r("idr.out$para$p"), ro.r("idr.out$para$rho")

    def get_p(p_val_1, p_val_2):
        return fit(p_val_1, p_val_2)[0]

    def get_rho(p_val_1, p_val_2):
        return fit(p_val_1, p_val_2)[1]
