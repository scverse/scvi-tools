import numpy as np

from scvi.data import synthetic_iid
from scvi.external import CellAssign


def test_cellassign(save_path):
    dataset = synthetic_iid(n_labels=5)
