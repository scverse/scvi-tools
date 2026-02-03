from __future__ import annotations

import logging
from collections import OrderedDict
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any


class FreezableEmbedding(nn.Embedding):
    """Embedding layer with partial freezing capability.

    This embedding layer allows freezing specific regions of the embedding
    matrix, which is useful for transfer learning and fine-tuning scenarios.

    Parameters
    ----------
    num_embeddings
        Number of embeddings (vocabulary size).
    embedding_dim
        Dimension of each embedding vector.
    n_freeze_x
        Number of embedding indices to freeze (from the beginning).
    n_freeze_y
        Number of embedding dimensions to freeze (from the beginning).
    **kwargs
        Additional arguments passed to nn.Embedding.

    Notes
    -----
    The freezing mechanism works by masking gradients during backpropagation.
    Frozen regions have their gradients set to zero, preventing parameter updates.

    Examples
    --------
    >>> import torch
    >>> # Create embedding with frozen first 5 indices and first 10 dimensions
    >>> emb = FreezableEmbedding(100, 50, n_freeze_x=5, n_freeze_y=10)
    >>> # Test forward pass
    >>> indices = torch.randint(0, 100, (10,))
    >>> output = emb(indices)
    >>> print(output.shape)  # torch.Size([10, 50])
    >>> # Test backward pass (frozen regions won't be updated)
    >>> loss = output.sum()
    >>> loss.backward()
    >>> print(emb.weight.grad[0, 0])  # Frozen region
    >>> print(emb.weight.grad[10, 20])  # Non-frozen region
    """

    def __init__(
        self, num_embeddings: int, embedding_dim: int, n_freeze_x: int = 0, n_freeze_y: int = 0, **kwargs: Any
    ) -> None:
        self._freeze_hook = None
        self.n_freeze_x = None
        self.n_freeze_y = None
        super().__init__(num_embeddings, embedding_dim, **kwargs)
        self.freeze(n_freeze_x, n_freeze_y)

    def freeze(self, n_freeze_x: int = 0, n_freeze_y: int = 0) -> None:
        """Set the freezing configuration for the embedding.

        Parameters
        ----------
        n_freeze_x
            Number of embedding indices to freeze (from the beginning).
        n_freeze_y
            Number of embedding dimensions to freeze (from the beginning).

        Notes
        -----
        This method registers a backward hook that masks gradients for
        the specified frozen regions. The hook is automatically removed
        and re-registered when this method is called.
        """
        self.n_freeze_x = n_freeze_x
        self.n_freeze_y = n_freeze_y

        if self._freeze_hook is not None:
            self._freeze_hook.remove()
            self._freeze_hook = None
        if self.n_freeze_x > 0 and self.n_freeze_y > 0:
            self._freeze_hook = self.weight.register_hook(self.partial_freeze_backward_hook)

    def partial_freeze_backward_hook(self, grad: torch.Tensor) -> torch.Tensor:
        """Backward hook to mask gradients for frozen regions.

        Parameters
        ----------
        grad
            Gradient tensor for the embedding weights.

        Returns
        -------
        torch.Tensor
            Masked gradient tensor.

        Notes
        -----
        This hook creates a mask that zeros out gradients for the frozen
        regions while preserving gradients for the non-frozen regions.
        """
        with torch.no_grad():
            assert self.n_freeze_x is not None and self.n_freeze_y is not None
            mask = F.pad(
                torch.zeros(self.n_freeze_x, self.n_freeze_y, device=grad.device),
                (0, self.embedding_dim - self.n_freeze_y, 0, self.num_embeddings - self.n_freeze_x),
                value=1.0,
            )
            return grad * mask

    def __repr__(self) -> str:
        """String representation of the embedding layer.

        Returns
        -------
        str
            String representation showing embedding dimensions and freeze status.
        """
        if self._freeze_hook is None:
            return f"Emb({self.num_embeddings}, {self.embedding_dim})"
        else:
            return f"Emb({self.num_embeddings}, {self.embedding_dim} | freeze: {self.n_freeze_x}, {self.n_freeze_y})"


class MultiEmbedding(nn.Module):
    """Multi-embedding layer that combines multiple embedding tables.

    This module combines multiple embedding layers with different vocabulary
    sizes and embedding dimensions into a single layer.

    Parameters
    ----------
    n_embedding_list
        List of vocabulary sizes for each embedding table.
    embedding_dim_list
        List of embedding dimensions for each embedding table.
    init_method
        Initialization method for embedding weights.
    normalization
        Normalization method to apply to concatenated embeddings.
    **kwargs
        Additional arguments passed to FreezableEmbedding.

    Notes
    -----
    The embeddings are concatenated along the last dimension. If normalization
    is specified, it's applied to the concatenated result.

    Examples
    --------
    >>> import torch
    >>> # Create multi-embedding with different vocab sizes and dimensions
    >>> multi_emb = MultiEmbedding(n_embedding_list=[100, 50, 200], embedding_dim_list=[32, 16, 64], normalization="l2")
    >>> # Test forward pass
    >>> indices = torch.randint(0, 100, (10, 3))  # 3 indices per sample
    >>> output = multi_emb(indices)
    >>> print(output.shape)  # torch.Size([10, 112])  # 32+16+64=112
    """

    def __init__(
        self,
        n_embedding_list: list[int],
        embedding_dim_list: list[int],
        init_method: str | Callable = "xavier_uniform",
        normalization: str | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__()
        assert len(n_embedding_list) == len(embedding_dim_list)

        self.emb_list = nn.ParameterList(
            [
                FreezableEmbedding(n_embedding, embedding_dim, **kwargs)
                for n_embedding, embedding_dim in zip(n_embedding_list, embedding_dim_list, strict=False)
            ]
        )
        assert normalization in [None, "l2"]
        self.normalization = normalization
        self.reset_parameters(init_method)

    def reset_parameters(self, init_method: str | Callable) -> None:
        """Reset parameters of all embedding layers.

        Parameters
        ----------
        init_method
            Initialization method to apply to each embedding layer.

        Notes
        -----
        Supported initialization methods:
        - "xavier_uniform": Xavier uniform initialization
        - "xavier_normal": Xavier normal initialization
        - "uniform": Uniform initialization in [-1, 1]
        - "normal": Normal initialization
        - "zero": Zero initialization
        - "one": One initialization
        - None: No initialization
        - callable: Custom initialization function
        """
        for emb in self.emb_list:
            if init_method is None:
                pass
            elif callable(init_method):
                init_method(emb.weight)
            elif init_method == "xavier_uniform":
                nn.init.xavier_uniform_(emb.weight)
            elif init_method == "xavier_normal":
                nn.init.xavier_normal_(emb.weight)
            elif init_method == "uniform":
                nn.init.uniform_(emb.weight, -1.0, 1.0)
            elif init_method == "normal":
                nn.init.normal_(emb.weight)
            elif init_method == "zero":
                nn.init.zeros_(emb.weight)
            elif init_method == "one":
                nn.init._no_grad_fill_(emb.weight, 1.0)
            else:
                raise NotImplementedError()

    def forward(self, index_list: torch.Tensor) -> torch.Tensor:
        """Forward pass through the multi-embedding layer.

        Parameters
        ----------
        index_list
            Tensor of indices with shape (..., n_embeddings).

        Returns
        -------
        torch.Tensor
            Concatenated embeddings with shape (..., total_embedding_dim).

        Notes
        -----
        The input should have the same number of indices as there are
        embedding tables. Each index is used to look up embeddings from
        the corresponding table, and the results are concatenated.
        """
        assert index_list.shape[-1] == len(self.emb_list)
        emb = torch.concat([emb(index_list[..., i]) for i, emb in enumerate(self.emb_list)], dim=-1)
        if self.normalization is None:
            return emb
        elif self.normalization == "l2":
            return F.normalize(emb, p=2, dim=1)
        else:
            raise NotImplementedError()

    @classmethod
    def from_pretrained(cls, feature_embedding_instance: Any) -> MultiEmbedding:
        """Create MultiEmbedding from a pretrained feature embedding.

        Parameters
        ----------
        feature_embedding_instance
            Pretrained feature embedding instance.

        Returns
        -------
        MultiEmbedding
            New MultiEmbedding instance with pretrained weights.

        Raises
        ------
        NotImplementedError
            This method is not yet implemented.
        """
        raise NotImplementedError()

    def load_weights_from_trained_module(self, other: MultiEmbedding, freeze_old: bool = False) -> None:
        """Load weights from another MultiEmbedding module.

        Parameters
        ----------
        other
            Source MultiEmbedding module to load weights from.
        freeze_old
            Whether to freeze the loaded weights.

        Notes
        -----
        This method transfers weights from a source module, handling cases
        where the target module has larger vocabulary sizes or embedding
        dimensions. The transferred weights are placed in the top-left
        corner of each embedding matrix.
        """
        assert len(self.emb_list) >= len(other.emb_list)
        if len(self.emb_list) > len(other.emb_list):
            logging.warning(f"Extending feature embedding {other} to {self} with more feature categories.")
        else:
            logging.info(f"Extending feature embedding {other} to {self}")
        for self_emb, other_emb in zip(self.emb_list, other.emb_list, strict=False):
            assert self_emb.num_embeddings >= other_emb.num_embeddings
            with torch.no_grad():
                extension_size = (
                    0,
                    self_emb.embedding_dim - other_emb.embedding_dim,
                    0,
                    self_emb.num_embeddings - other_emb.num_embeddings,
                )
                transfer_mask = F.pad(
                    torch.zeros_like(other_emb.weight.data, device=other_emb.weight.data.device),
                    extension_size,
                    value=1.0,
                )
                extended_other_emb = F.pad(other_emb.weight.data, extension_size, value=0.0)
                self_emb.weight.data = extended_other_emb + transfer_mask * self_emb.weight.data
            if freeze_old:
                self_emb.freeze(other_emb.num_embeddings, other_emb.embedding_dim)

    def freeze_top_embs(self, n_freeze_list: list[int]) -> None:
        """Freeze the top embeddings of each embedding table.

        Parameters
        ----------
        n_freeze_list
            Number of top embeddings to freeze for each table.

        Notes
        -----
        This method freezes the first n embeddings from each embedding
        table, which is useful for transfer learning scenarios.
        """
        for emb, n_freeze in zip(self.emb_list, n_freeze_list, strict=False):
            # If that specific category has no change in size, skip.
            if not emb.weight.requires_grad:
                continue
            emb.freeze(n_freeze, emb.embedding_dim)

    @property
    def num_embeddings(self) -> list[int]:
        """Get list of vocabulary sizes.

        Returns
        -------
        list[int]
            List of vocabulary sizes for each embedding table.
        """
        return [emb.num_embeddings for emb in self.emb_list]

    @property
    def embedding_dim(self) -> int:
        """Get total embedding dimension.

        Returns
        -------
        int
            Sum of all embedding dimensions.
        """
        return sum(emb.embedding_dim for emb in self.emb_list)

    def __repr__(self) -> str:
        """String representation of the multi-embedding layer.

        Returns
        -------
        str
            String representation showing all embedding tables and normalization.
        """
        repr_text = "cat(" + ", ".join(repr(emb) for emb in self.emb_list) + ")"
        if self.normalization:
            repr_text = f"{self.normalization}({repr_text})"
        return repr_text


class FeatureEmbedding(nn.Module):
    """Feature embedding layer for categorical features.

    This module creates embeddings for categorical features based on
    vocabulary lists. It handles the mapping from categorical values
    to embedding indices automatically.

    Parameters
    ----------
    vocab_list
        List of vocabulary lists, one for each categorical feature.
    embedding_dims
        List of embedding dimensions for each feature.
    **kwargs
        Additional arguments passed to MultiEmbedding.

    Notes
    -----
    This module automatically creates vocabulary mappings and handles
    the conversion from categorical values to indices. It also provides
    caching for efficient repeated lookups.

    Examples
    --------
    >>> import numpy as np
    >>> # Create feature embedding for categorical features
    >>> vocab_list = [["A", "B", "C"], ["X", "Y"], ["1", "2", "3", "4"]]
    >>> embedding_dims = [16, 8, 32]
    >>> feature_emb = FeatureEmbedding(vocab_list, embedding_dims)
    >>> # Test with categorical data
    >>> sentences = np.array([["A", "X", "1"], ["B", "Y", "2"]])
    >>> output = feature_emb(sentences)
    >>> print(output.shape)  # torch.Size([2, 56])  # 16+8+32=56
    """

    def __init__(self, vocab_list: list[list[str]], embedding_dims: list[int], **kwargs: Any) -> None:
        super().__init__()
        assert len(vocab_list) == len(embedding_dims)
        self.device_container = nn.Parameter(torch.tensor([]))

        self.vocab_list = vocab_list
        n_vocab_list = [len(vocab) for vocab in self.vocab_list]
        self.multi_emb = self.define_embeddings(n_vocab_list, embedding_dims, **kwargs)

        self.__index_cache: dict[str, torch.Tensor] = {}
        self.__vocab_map_list_cache: list[dict[str, int]] | None = None

    def reset_cache(self) -> None:
        """Reset the internal caches.

        Notes
        -----
        This method clears both the index cache and vocabulary mapping cache.
        It should be called when the vocabulary changes.
        """
        self.__index_cache = {}
        self.__vocab_map_list_cache = None

    @property
    def vocab_map_list(self) -> list[dict[str, int]]:
        """Get list of vocabulary mappings.

        Returns
        -------
        list[dict]
            List of dictionaries mapping categorical values to indices.

        Notes
        -----
        This property is cached for efficiency. The cache is automatically
        invalidated when reset_cache() is called.
        """
        if self.__vocab_map_list_cache is None:
            n_vocab_list = [len(vocab) for vocab in self.vocab_list]
            self.__vocab_map_list_cache = [
                dict(zip(vocab, range(n_vocab), strict=False))
                for vocab, n_vocab in zip(self.vocab_list, n_vocab_list, strict=False)
            ]
        return self.__vocab_map_list_cache

    @staticmethod
    def define_embeddings(n_embedding_list: list[int], embedding_dims: list[int], **kwargs: Any) -> MultiEmbedding:
        """Define the underlying embedding layers.

        Parameters
        ----------
        n_embedding_list
            List of vocabulary sizes.
        embedding_dims
            List of embedding dimensions.
        **kwargs
            Additional arguments for MultiEmbedding.

        Returns
        -------
        MultiEmbedding
            Multi-embedding layer.
        """
        return MultiEmbedding(n_embedding_list, embedding_dims, **kwargs)

    def reset_parameters(self, init_method: str | Callable) -> None:
        """Reset parameters of the embedding layers.

        Parameters
        ----------
        init_method : str or callable
            Initialization method to apply.
        """
        self.multi_emb.reset_parameters(init_method)

    def _get_index_from_sentences(self, index_sentences: np.ndarray) -> torch.Tensor:
        """Convert categorical sentences to embedding indices.

        Parameters
        ----------
        index_sentences
            Array of categorical values with shape (..., n_features).

        Returns
        -------
        torch.Tensor
            Tensor of embedding indices with shape (..., n_features).

        Notes
        -----
        This method uses the vocabulary mappings to convert categorical
        values to their corresponding embedding indices.
        """
        assert index_sentences.shape[-1] == len(self.vocab_map_list)
        mapping_list = map(lambda mapping: np.vectorize(lambda key: mapping[key]), self.vocab_map_list)

        indices = torch.concat(
            [torch.from_numpy(mapping(index_sentences[..., [i]])) for i, mapping in enumerate(mapping_list)], dim=-1
        )
        return indices.to(self.device_container.device)

    def forward(self, index_sentences: np.ndarray, index_cache_key: str | None = None) -> torch.Tensor:
        """Forward pass through the feature embedding layer.

        Parameters
        ----------
        index_sentences
            Array of categorical values with shape (..., n_features).
        index_cache_key
            Cache key for storing computed indices.

        Returns
        -------
        torch.Tensor
            Concatenated embeddings with shape (..., total_embedding_dim).

        Notes
        -----
        If a cache key is provided, the computed indices are stored for
        future use, improving efficiency for repeated lookups.
        """
        if index_cache_key is not None and index_cache_key in self.__index_cache:
            return self.multi_emb(self.__index_cache[index_cache_key])
        indices = self._get_index_from_sentences(index_sentences)
        if index_cache_key is not None:
            self.__index_cache[index_cache_key] = indices
        return self.multi_emb(indices)

    def get_extra_state(self) -> dict[str, Any]:
        """Get extra state for serialization.

        Returns
        -------
        dict
            Extra state information including vocabulary lists.
        """
        return {"vocab_list": self.vocab_list}

    def set_extra_state(self, state: dict[str, Any]) -> None:
        """Set extra state from serialization.

        Parameters
        ----------
        state
            Extra state information.
        """
        self.vocab_list = state["vocab_list"]
        self.__vocab_map_list_cache = None

    @classmethod
    def from_numpy_array(
        cls, sentences_array: np.ndarray, embedding_dims: list[int], **kwargs: Any
    ) -> FeatureEmbedding:
        """Create FeatureEmbedding from a numpy array.

        Parameters
        ----------
        sentences_array
            Array of categorical sentences with shape (n_samples, n_features).
        embedding_dims
            Embedding dimensions for each feature.
        **kwargs
            Additional arguments for FeatureEmbedding.

        Returns
        -------
        FeatureEmbedding
            New FeatureEmbedding instance.

        Notes
        -----
        This method automatically extracts unique values for each feature
        to create the vocabulary lists.
        """
        word_list = sentences_array.transpose().tolist()
        vocab_list = [list(OrderedDict.fromkeys(words)) for words in word_list]
        return cls(vocab_list, embedding_dims, **kwargs)

    @classmethod
    def from_pandas_dataframe(
        cls, sentences_df: pd.DataFrame, embedding_dims: list[int], **kwargs: Any
    ) -> FeatureEmbedding:
        """Create FeatureEmbedding from a pandas DataFrame.

        Parameters
        ----------
        sentences_df
            DataFrame of categorical sentences.
        embedding_dims
            Embedding dimensions for each feature.
        **kwargs
            Additional arguments for FeatureEmbedding.

        Returns
        -------
        FeatureEmbedding
            New FeatureEmbedding instance.
        """
        if isinstance(embedding_dims, pd.DataFrame):
            assert sentences_df.columns == embedding_dims.columns
            assert len(embedding_dims) == 1
            embedding_dims = embedding_dims.loc[0].to_list()
        return cls.from_numpy_array(sentences_df.values, embedding_dims, **kwargs)

    @classmethod
    def from_pretrained(cls, feature_embedding_instance: FeatureEmbedding) -> FeatureEmbedding:
        """Create FeatureEmbedding from a pretrained instance.

        Parameters
        ----------
        feature_embedding_instance
            Pretrained FeatureEmbedding instance.

        Returns
        -------
        FeatureEmbedding
            New FeatureEmbedding instance with pretrained weights.

        Raises
        ------
        NotImplementedError
            This method is not yet implemented.
        """
        raise NotImplementedError()

    def load_weights_from_trained_module(self, other: FeatureEmbedding, freeze_old: bool = False) -> None:
        """Load weights from another FeatureEmbedding module.

        Parameters
        ----------
        other
            Source FeatureEmbedding module to load weights from.
        freeze_old
            Whether to freeze the loaded weights.

        Notes
        -----
        This method transfers weights from the source module's multi_emb
        to this module's multi_emb.
        """
        assert isinstance(other, self.__class__)
        assert len(self.vocab_list) >= len(other.vocab_list)

        for self_vocab, other_vocab in zip(self.vocab_list, other.vocab_list, strict=False):
            assert len(self_vocab) >= len(other_vocab)
            for self_word, other_word in zip(self_vocab, other_vocab, strict=False):
                assert self_word == other_word

        self.multi_emb.load_weights_from_trained_module(other.multi_emb, freeze_old=freeze_old)

    @property
    def embedding_dim(self) -> int:
        """Get total embedding dimension.

        Returns
        -------
        int
            Sum of all embedding dimensions.
        """
        return self.multi_emb.embedding_dim

    def __repr__(self) -> str:
        """String representation of the feature embedding layer.

        Returns
        -------
        str
            String representation showing vocabulary sizes and embedding dimensions.
        """
        return f"str2index -> {repr(self.multi_emb)}"
