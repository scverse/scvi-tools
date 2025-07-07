import logging

import anndata
import numpy as np
from scipy import sparse

logger = logging.getLogger(__name__)


class ArrayFakerSparse(sparse.csr_matrix):
    """Create a wrapper around a sparse matrix layer in an AnnCollection.

    This wrapper ensures compatibility with the typical AnnData API when
    accessing or manipulating sparse matrix layers.
    """

    def __init__(self, collection, key):
        self.collection = collection
        self.key = key
        # minibatch size for batchwise operations
        self.batch_size = 256

    def __getitem__(self, idx):
        """Second layer of `__getitem__` that handles returning an index into a specific layer."""
        return self.collection[idx].layers[self.key]

    def __repr__(self):
        return f"ArrayFakerSparse(collection={self.collection}, key={self.key})"

    @property
    def data(self):
        """Implement a reference into a subset of the non-zero values.

        Notes
        -----
        This is used in `scvi-tools` to check that values are integers. We
        just return a subsample for now so that existing code paths are
        supported.

        """
        return self.collection[: np.min([64, len(self.collection)])].layers[self.key].data

    def getformat(self):
        """Extract the relevant format (CSR, CSC, COO)"""
        return self.collection[: np.min([4, len(self.collection)])].layers[self.key].getformat()


class LayerFaker:
    """Custom access for .layers inside an AnnCollection.

    This class provides support for accessing layers when working with
    AnnCollection-backed datasets through a CollectionAdapter.
    """

    def __init__(self, collection):
        self.collection = collection

    def keys(self):
        k = self.collection[1 : np.min([4, len(self.collection)])].layers.keys()
        return k

    def __getitem__(self, key):
        """First layer of get item that handles returning an index into a specific layer"""
        # determine the array's type
        # we first need to draw an AnnCollectionView to get an object that has
        # access to .layers, which is not provided in AnnCollection
        sample = self.collection[: np.min([4, len(self.collection)])]
        layer = sample.layers[key]

        if isinstance(layer, sparse.spmatrix):
            return ArrayFakerSparse(self.collection, key)
        elif isinstance(layer, np.ndarray):
            # TODO: Not yet implemented
            raise TypeError("ArrayFakerNDarray is not yet implemented")
        else:
            raise TypeError(f"Unknown array type {type(layer)}")


class CollectionAdapter:
    """Allow an AnnCollection to pretend to be an AnnData in SCVI-Tools"""

    # scvi-tools stores registry information here
    uns = {}
    # necessary for an scvi-tools check
    isbacked = True
    # necessary to pass an scvi-tools check that rejects views
    is_view = False

    SPECIAL_ATTRS = [
        "uns",
        "is_view",
        "isbacked",
        "layers",
    ]

    def __init__(self, collection):
        self.collection = collection

    def __getattr__(self, name):
        if name in self.SPECIAL_ATTRS:
            return getattr(self, name)
        else:
            return getattr(self.collection, name)

    def __repr__(self):
        return "Adapter for:\n" + repr(self.collection)

    @property
    def layers(self):
        return LayerFaker(self.collection)

    def __getitem__(self, idx):
        return self.collection[idx]

    def __len__(self):
        return len(self.collection)

    @property
    def __class__(self):
        return anndata.AnnData
