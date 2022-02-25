from typing import Optional, Union

import numpy as np
import pandas as pd
import torch
from scipy.sparse import spmatrix
from torch import Tensor, isin

from scvi._compat import Literal


def mde(
    data: Union[np.ndarray, pd.DataFrame, spmatrix, Tensor],
    device: Optional[Literal["cpu", "cuda"]] = None,
    **kwargs,
) -> np.ndarray:
    """
    Util to run :func:`pymde.preserve_neighbors` for visualization of scvi-tools embeddings.

    Parameters
    ----------
    data
        The data of shape (n_obs, k), where k is typically defined by one of the models
        in scvi-tools that produces an embedding (e.g., :class:`~scvi.model.SCVI`.)
    device
        Whether to run on cpu or gpu ("cuda"). If None, tries to run on gpu if available.
    kwargs
        Keyword args to :func:`pymde.preserve_neighbors`


    Notes
    -----
    This function is included in scvi-tools to provide an alternative to UMAP/TSNE that is GPU-
    accelerated. The appropriateness of use of visualization of high-dimensional spaces in single-
    cell omics remains an open research questions. See:

    Chari, Tara, Joeyta Banerjee, and Lior Pachter. "The specious art of single-cell genomics." bioRxiv (2021).

    If you use this function in your research please cite:

    Agrawal, Akshay, Alnur Ali, and Stephen Boyd. "Minimum-distortion embedding." arXiv preprint arXiv:2103.02559 (2021).

    Returns
    -------
    The pymde embedding, defaults to two dimensions.
    """
    try:
        import pymde
    except ImportError:
        raise ImportError("Please install pymde package via `pip install pymde`")

    if isinstance(data, pd.DataFrame):
        data = data.values

    device = "cpu" if not torch.cuda.is_available() else "cuda"

    _kwargs = dict(
        embedding_dim=2,
        constraint=pymde.Standardized(),
        repulsive_fraction=0.7,
        verbose=False,
        device=device,
    )
    _kwargs.update(kwargs)

    emb = pymde.preserve_neighbors(data, **_kwargs).embed(verbose=_kwargs["verbose"])

    if isinstance(emb, torch.Tensor):
        emb = emb.cpu().numpy()

    return emb
