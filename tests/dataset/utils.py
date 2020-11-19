import torch

from anndata import AnnData

import scvi

use_gpu = torch.cuda.is_available()


def unsupervised_training_one_epoch(adata: AnnData):
    m = scvi.model.SCVI(adata)
    m.train(1, train_size=0.4, use_gpu=use_gpu)
