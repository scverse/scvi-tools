import pytest

import scvi
from scvi.data.fields import LabelsWithUnlabeledObsField
from scvi.model import SCANVI


def test_scanvi_unlabeled_category(n_labels: int = 3):
    adata = scvi.data.synthetic_iid(n_labels=n_labels)

    # unlabeled category passed in
    SCANVI.setup_anndata(adata, labels_key="labels", unlabeled_category="label_0")
    model = SCANVI(adata)
    model.train(max_epochs=1)

    # unlabeled category not passed in
    SCANVI.setup_anndata(adata, labels_key="labels")
    model = SCANVI(adata)
    model.train(max_epochs=1)

    # unlabeled category passed in, but not present in labels
    # with pytest.raises(ValueError):
    #     SCANVI.setup_anndata(
    #         adata, labels_key="labels", unlabeled_category="invalid_label"
    #     )
    # TODO: currently does not raise error but should

    # unlabeled category not passed in, but labels contains reserved category
    adata.obs["labels"] = LabelsWithUnlabeledObsField.UNLABELED_CATEGORY_KEY
    with pytest.raises(ValueError):
        SCANVI.setup_anndata(adata, labels_key="labels")
