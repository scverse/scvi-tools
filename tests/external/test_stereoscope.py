import scvi
import os

import numpy as np
from scvi.external.stereoscope import RNAStereoscope, SpatialStereoscope
from scvi.data import register_tensor_from_anndata

dataset = scvi.data.pbmc_dataset(
            save_path="tests/data/10X",
            remove_extracted_data=True,
            run_setup_anndata=True,
        )


model = RNAStereoscope(dataset)
model.train(n_epochs=100, frequency=1, vae_task_kwargs={"lr":0.01}, train_size=1.0)
params = model.get_params()
print(model.history)


dataset.obs["indices"] = np.arange(dataset.n_obs)
register_tensor_from_anndata(dataset, "ind_x", "obs", "indices")

model = SpatialStereoscope(dataset, params)
model.train(n_epochs=100, frequency=1, vae_task_kwargs={"lr":0.01}, train_size=1.0)
print(model.get_proportions())