from typing import Optional

import torch
from anndata import AnnData

from scvi.model import SCVI

use_gpu = torch.cuda.is_available()


def unsupervised_training_one_epoch(
    adata: AnnData,
    run_setup_anndata: bool = True,
    batch_key: Optional[str] = None,
    labels_key: Optional[str] = None,
):
    if run_setup_anndata:
        SCVI.setup_anndata(adata, batch_key=batch_key, labels_key=labels_key)
    m = SCVI(adata)
    m.train(1, train_size=0.4, use_gpu=use_gpu)
