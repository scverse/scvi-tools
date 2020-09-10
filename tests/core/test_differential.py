from scvi.dataset import synthetic_iid
from scvi.models import SCVI
from scvi.models.differential import DifferentialComputation
from functools import partial
import numpy as np


def test_differential_computation(save_path):

    n_latent = 5
    adata = synthetic_iid()
    model = SCVI(adata, n_latent=n_latent)
    model.train(1)

    model_fn = partial(model.get_normalized_expression, return_numpy=True)
    dc = DifferentialComputation(model_fn, adata)

    cell_idx1 = np.asarray(adata.obs.labels == "undefined_1")
    cell_idx2 = ~cell_idx1

    dc.get_bayes_factors(cell_idx1, cell_idx2, mode="vanilla", use_permutation=True)
    dc.get_bayes_factors(cell_idx1, cell_idx2, mode="change", use_permutation=False)
    dc.get_bayes_factors(cell_idx1, cell_idx2, mode="change", cred_interval_lvls=[0.75])

    delta = 0.5

    def change_fn_test(x, y):
        return x - y

    def m1_domain_fn_test(samples):
        return np.abs(samples) >= delta

    dc.get_bayes_factors(
        cell_idx1,
        cell_idx2,
        mode="change",
        m1_domain_fn=m1_domain_fn_test,
        change_fn=change_fn_test,
    )
