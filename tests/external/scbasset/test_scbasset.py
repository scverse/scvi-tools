from pathlib import Path

import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.external import SCBASSET

_DNA_CODE_KEY = "code"


def _get_adata(sparse_format: str | None = None):
    dataset1 = synthetic_iid(batch_size=100, sparse_format=sparse_format).transpose()
    dataset1.X = (dataset1.X > 0).astype(float)
    dataset1.obsm[_DNA_CODE_KEY] = np.random.randint(0, 3, size=(dataset1.n_obs, 1344))
    return dataset1


def test_scbasset():
    adata = _get_adata()
    SCBASSET.setup_anndata(
        adata,
        dna_code_key=_DNA_CODE_KEY,
    )
    model = SCBASSET(adata)
    model.train(max_epochs=2, early_stopping=True)
    model.get_latent_representation()


def test_scbasset_batch():
    adata = _get_adata()
    SCBASSET.setup_anndata(
        adata,
        dna_code_key=_DNA_CODE_KEY,
        batch_key="batch",
    )
    model = SCBASSET(adata)
    model.train(max_epochs=1)
    model.get_latent_representation()
    assert hasattr(model.module, "batch_ids")


@pytest.mark.internet
def test_scbasset_motif_download(save_path):
    # get a temporary directory name
    SCBASSET._download_motifs(genome="human", motif_dir=save_path)

    assert Path(save_path, "shuffled_peaks.fasta").exists()
    assert Path(save_path, "shuffled_peaks_motifs", "MYOD1.fasta").exists()
    return


@pytest.mark.internet
def test_scbasset_tf(save_path):
    adata = _get_adata()
    SCBASSET.setup_anndata(
        adata,
        dna_code_key=_DNA_CODE_KEY,
    )
    model = SCBASSET(adata)
    # model.train(max_epochs=2, early_stopping=True)
    model.is_trained = True

    motif_seqs, bg_seqs = model._get_motif_library(tf="MYOD1", motif_dir=save_path, genome="human")

    model.get_tf_activity(tf="MYOD1", motif_dir=save_path, genome="human", lib_size_norm=True)
    model.get_tf_activity(tf="MYOD1", motif_dir=save_path, genome="human", lib_size_norm=False)
    return
