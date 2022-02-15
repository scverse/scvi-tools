import numpy as np

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

    st_model = SpatialStereoscope.from_rna_model(
        dataset, sc_model, prior_weight="minibatch"
    )
    st_model.train(max_epochs=1)
    st_model.get_proportions()
    # test save/load
    st_model.save(save_path, overwrite=True, save_anndata=True)
    st_model = SpatialStereoscope.load(save_path)
    st_model.get_proportions()

    # try imputation code
    y = np.array(50 * ["label_0"])
    st_model.get_scale_for_ct(y)
