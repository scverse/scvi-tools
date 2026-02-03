from __future__ import annotations

from collections import OrderedDict
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F

from scvi.external.drvi.nn_modules.embedding import FeatureEmbedding

if TYPE_CHECKING:
    from typing import Any


class MultiOneHotEncoding(nn.Module):
    """Multi-one-hot encoding layer for categorical features.

    This module creates one-hot encodings for multiple categorical features
    and concatenates them. It's similar to MultiEmbedding but uses one-hot
    encoding instead of learned embeddings.

    Parameters
    ----------
    n_embedding_list
        List of vocabulary sizes for each categorical feature.
    **kwargs
        Additional arguments (unused, kept for compatibility).

    Notes
    -----
    This module creates one-hot encodings for each categorical feature and
    concatenates them along the last dimension. The total output dimension
    is the sum of all vocabulary sizes.

    Examples
    --------
    >>> import torch
    >>> # Create multi-one-hot encoding
    >>> encoding = MultiOneHotEncoding([3, 2, 4])
    >>> # Test with indices
    >>> indices = torch.tensor([[0, 1, 2], [2, 0, 1]])
    >>> output = encoding(indices)
    >>> print(output.shape)  # torch.Size([2, 9])  # 3+2+4=9
    >>> # Check one-hot encoding
    >>> print(output[0])  # [1,0,0, 0,1, 0,0,1,0] for indices [0,1,2]
    """

    def __init__(self, n_embedding_list: list[int], **kwargs: Any) -> None:
        super().__init__()
        self.n_embedding_list = n_embedding_list
        self.device_container = nn.Parameter(torch.tensor([]))

    def forward(self, index_list: torch.Tensor) -> torch.Tensor:
        """Forward pass through the multi-one-hot encoding layer.

        Parameters
        ----------
        index_list
            Tensor of indices with shape (..., n_features).

        Returns
        -------
        torch.Tensor
            Concatenated one-hot encodings with shape (..., total_vocab_size).

        Notes
        -----
        Each index is converted to a one-hot encoding, and all encodings
        are concatenated along the last dimension. The output dimension
        is the sum of all vocabulary sizes.
        """
        assert index_list.shape[-1] == len(self.n_embedding_list)
        result = torch.concat(
            [
                F.one_hot(index_list[..., i], num_classes=n_embedding)
                for i, n_embedding in enumerate(self.n_embedding_list)
            ],
            dim=-1,
        )
        return result.to(self.device_container.device)

    def get_extra_state(self) -> dict[str, Any]:
        """Get extra state for serialization.

        Returns
        -------
        dict
            Extra state information including vocabulary sizes.
        """
        return {"n_embedding_list": self.n_embedding_list}

    def set_extra_state(self, state: dict[str, Any]) -> None:
        """Set extra state from serialization.

        Parameters
        ----------
        state
            Extra state information.
        """
        self.n_embedding_list = state["n_embedding_list"]

    @property
    def embedding_dim(self) -> int:
        """Get total embedding dimension.

        Returns
        -------
        int
            Sum of all vocabulary sizes.

        Notes
        -----
        For one-hot encoding, the embedding dimension equals the vocabulary size
        since each category is represented by a single dimension.
        """
        return sum(self.n_embedding_list)


class FeatureOneHotEncoding(FeatureEmbedding):
    """Feature one-hot encoding layer for categorical features.

    This module creates one-hot encodings for categorical features based on
    vocabulary lists. It's similar to FeatureEmbedding but uses one-hot
    encoding instead of learned embeddings.

    Parameters
    ----------
    vocab_list
        List of vocabulary lists, one for each categorical feature.
    **kwargs
        Additional arguments passed to FeatureEmbedding.

    Notes
    -----
    This module automatically creates vocabulary mappings and handles
    the conversion from categorical values to one-hot encodings. It
    inherits most functionality from FeatureEmbedding but uses one-hot
    encoding instead of learned embeddings.

    Examples
    --------
    >>> import numpy as np
    >>> # Create feature one-hot encoding
    >>> vocab_list = [["A", "B", "C"], ["X", "Y"], ["1", "2", "3", "4"]]
    >>> encoding = FeatureOneHotEncoding(vocab_list)
    >>> # Test with categorical data
    >>> sentences = np.array([["A", "X", "1"], ["B", "Y", "2"]])
    >>> output = encoding(sentences)
    >>> print(output.shape)  # torch.Size([2, 9])  # 3+2+4=9
    >>> # Check one-hot encoding for first sample
    >>> print(output[0])  # One-hot encoding for ["A", "X", "1"]
    """

    def __init__(self, vocab_list: list[list[str]], **kwargs: Any) -> None:
        n_vocab_list = [len(vocab) for vocab in vocab_list]
        super().__init__(vocab_list, n_vocab_list, **kwargs)

    @staticmethod
    def define_embeddings(n_embedding_list: list[int], embedding_dims: list[int], **kwargs: Any) -> MultiOneHotEncoding:
        """Define the underlying one-hot encoding layers.

        Parameters
        ----------
        n_embedding_list
            List of vocabulary sizes.
        embedding_dims
            List of embedding dimensions (must equal vocabulary sizes).
        **kwargs
            Additional arguments for MultiOneHotEncoding.

        Returns
        -------
        MultiOneHotEncoding
            Multi-one-hot encoding layer.

        Notes
        -----
        For one-hot encoding, the embedding dimension must equal the
        vocabulary size since each category is represented by a single
        dimension in the one-hot vector.
        """
        for n_key, dim in zip(n_embedding_list, embedding_dims, strict=False):
            assert n_key == dim
        return MultiOneHotEncoding(n_embedding_list, **kwargs)

    @classmethod
    def from_numpy_array(cls, sentences_array: np.ndarray, **kwargs: Any) -> FeatureOneHotEncoding:
        """Create FeatureOneHotEncoding from a numpy array.

        Parameters
        ----------
        sentences_array
            Array of categorical sentences with shape (n_samples, n_features).
        **kwargs
            Additional arguments for FeatureOneHotEncoding.

        Returns
        -------
        FeatureOneHotEncoding
            New FeatureOneHotEncoding instance.

        Notes
        -----
        This method automatically extracts unique values for each feature
        to create the vocabulary lists. The embedding dimensions are set
        to match the vocabulary sizes for one-hot encoding.
        """
        word_list = sentences_array.transpose().tolist()
        vocab_list = [list(OrderedDict.fromkeys(words)) for words in word_list]
        return cls(vocab_list, **kwargs)

    @classmethod
    def from_pandas_dataframe(cls, sentences_df: pd.DataFrame, **kwargs: Any) -> FeatureOneHotEncoding:
        """Create FeatureOneHotEncoding from a pandas DataFrame.

        Parameters
        ----------
        sentences_df : pandas.DataFrame
            DataFrame of categorical sentences.
        **kwargs
            Additional arguments for FeatureOneHotEncoding.

        Returns
        -------
        FeatureOneHotEncoding
            New FeatureOneHotEncoding instance.
        """
        return cls.from_numpy_array(sentences_df.values, **kwargs)
