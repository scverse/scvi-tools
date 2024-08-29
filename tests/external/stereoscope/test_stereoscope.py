from time import time

import numpy as np
import torch

from scvi.data import synthetic_iid
from scvi.external import RNAStereoscope, SpatialStereoscope


def test_stereoscope(save_path):
    dataset = synthetic_iid(
        n_labels=5,
    )
    RNAStereoscope.setup_anndata(
        dataset,
        labels_key="labels",
    )

    # train with no proportions
    sc_model = RNAStereoscope(dataset)
    sc_model.train(max_epochs=1)

    # train again with proportions
    sc_model = RNAStereoscope(dataset, ct_weights=np.ones((5,)))
    sc_model.train(max_epochs=1)
    # test save/load
    sc_model.save(save_path, overwrite=True, save_anndata=True)
    sc_model = RNAStereoscope.load(save_path)
    dataset = synthetic_iid(
        n_labels=5,
    )
    SpatialStereoscope.setup_anndata(
        dataset,
    )
    st_model = SpatialStereoscope.from_rna_model(dataset, sc_model, prior_weight="minibatch")
    st_model.train(max_epochs=1)
    st_model.get_proportions()
    # test save/load
    st_model.save(save_path, overwrite=True, save_anndata=True)
    st_model = SpatialStereoscope.load(save_path)
    st_model.get_proportions()

    # try imputation code
    y = np.array(50 * ["label_0"])
    st_model.get_scale_for_ct(y)


def test_cpu_gpu_stereoscope():
    if torch.cuda.is_available():
        adata = synthetic_iid(10000, 500)

        RNAStereoscope.setup_anndata(adata, labels_key="labels")

        m = RNAStereoscope(adata)
        training_start_time = time()
        m.train(
            accelerator="cpu",
            batch_size=5000,
            max_epochs=100,
            train_size=0.9,
            plan_kwargs={"n_epochs_kl_warmup": 100, "compile": False},
            datasplitter_kwargs={"drop_last": True},
        )
        print(f"CPU Training finished, took {time() - training_start_time:.2f}s")

        # run the exact same thing on GPU:
        m2 = RNAStereoscope(adata)
        training_start_time2 = time()
        m2.train(
            accelerator="cuda",
            batch_size=5000,
            max_epochs=100,
            train_size=0.9,
            plan_kwargs={"n_epochs_kl_warmup": 100, "compile": True},
            datasplitter_kwargs={"drop_last": True},
        )
        print(f"Compile + GPU Training finished, took {time() - training_start_time2:.2f}s")
