import numpy as np

from scvi.data import register_tensor_from_anndata, synthetic_iid
from scvi.external import RNAStereoscope, SpatialStereoscope


def test_stereoscope(save_path):
    dataset = synthetic_iid()
    sc_model = RNAStereoscope(dataset)
    sc_model.train(
        n_epochs=1, frequency=1, vae_task_kwargs={"lr": 0.01}, train_size=1.0
    )
    sc_params = sc_model.get_params()

    dataset.obs["indices"] = np.arange(dataset.n_obs)
    register_tensor_from_anndata(dataset, "ind_x", "obs", "indices")

    st_model = SpatialStereoscope(dataset, sc_params)
    st_model.train(
        n_epochs=1, frequency=1, vae_task_kwargs={"lr": 0.01}, train_size=1.0
    )
    st_model.get_proportions()
