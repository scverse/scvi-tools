from __future__ import annotations

import abc
import gc
import logging
import os
import threading
from collections import deque
from collections.abc import Iterator, Sequence
from concurrent import futures
from concurrent.futures import Future
from contextlib import contextmanager
from datetime import timedelta
from math import ceil
from time import time
from typing import Any, TypeVar

import numpy as np
import numpy.typing as npt
import pandas as pd
import psutil
import scipy
import tiledbsoma as soma
import torch
import torchdata.datapipes.iter as pipes
from attr import define
from lightning.pytorch import LightningDataModule
from numpy.random import Generator
from scipy import sparse
from sklearn.preprocessing import LabelEncoder
from torch import Tensor
from torch import distributed as dist
from torch.utils.data import DataLoader
from torch.utils.data.dataset import Dataset

pytorch_logger = logging.getLogger("cellxgene_census.experimental.pytorch")

# TODO: Rename to reflect the correct order of the Tensors within the tuple: (X, obs)
ObsAndXDatum = tuple[Tensor, Tensor]
"""Return type of ``ExperimentDataPipe`` that pairs a Tensor of ``obs`` row(s) with a Tensor of
``X`` matrix row(s).The Tensors are rank 1 if ``batch_size`` is 1,
otherwise the Tensors are rank 2."""

util_logger = logging.getLogger("cellxgene_census.experimental.util")

_T = TypeVar("_T")


DEFAULT_TILEDB_CONFIGURATION: dict[str, Any] = {
    # https://docs.tiledb.com/main/how-to/configuration#configuration-parameters
    "py.init_buffer_bytes": 1 * 1024**3,
    "soma.init_buffer_bytes": 1 * 1024**3,
    # S3 requests should not be signed, since we want to allow anonymous access
    "vfs.s3.no_sign_request": "true",
    "vfs.s3.region": "us-west-2",
}


def get_default_soma_context(
    tiledb_config: dict[str, Any] | None = None,
) -> soma.options.SOMATileDBContext:
    """Return a :class:`tiledbsoma.SOMATileDBContext` with sensible defaults that can be further

    customized by the user. The customized context can then be passed to
    :func:`cellxgene_census.open_soma` with the ``context`` argument or to
    :meth:`somacore.SOMAObject.open` with the ``context`` argument, such as
    :meth:`tiledbsoma.Experiment.open`. Use the :meth:`tiledbsoma.SOMATileDBContext.replace`
    method on the returned object to customize its settings further.

    Args:
        tiledb_config:
            A dictionary of TileDB configuration parameters. If specified, the parameters will
            override the defaults. If not specified, the default configuration will be returned.

    Returns
    -------
        A :class:`tiledbsoma.SOMATileDBContext` object with sensible defaults.

    Examples
    --------
        To reduce the amount of memory used by TileDB-SOMA I/O operations:

        .. highlight:: python
        .. code-block:: python

            ctx = cellxgene_census.get_default_soma_context(
                tiledb_config={
                    "py.init_buffer_bytes": 128 * 1024**2,
                    "soma.init_buffer_bytes": 128 * 1024**2,
                }
            )
            c = census.open_soma(uri="s3://my-private-bucket/census/soma", context=ctx)

        To access a copy of the Census located in a private bucket that is located in a different
        S3 region, use:

        .. highlight:: python
        .. code-block:: python

            ctx = cellxgene_census.get_default_soma_context(
                tiledb_config={"vfs.s3.no_sign_request": "false", "vfs.s3.region": "us-east-1"}
            )
            c = census.open_soma(uri="s3://my-private-bucket/census/soma", context=ctx)

    Lifecycle:
        experimental
    """
    tiledb_config = dict(DEFAULT_TILEDB_CONFIGURATION, **(tiledb_config or {}))
    return soma.options.SOMATileDBContext().replace(tiledb_config=tiledb_config)


class _EagerIterator(Iterator[_T]):
    def __init__(
        self,
        iterator: Iterator[_T],
        pool: futures.Executor | None = None,
    ):
        super().__init__()
        self.iterator = iterator
        self._pool = pool or futures.ThreadPoolExecutor()
        self._own_pool = pool is None
        self._future: Future[_T] | None = None
        self._begin_next()

    def _begin_next(self) -> None:
        self._future = self._pool.submit(self.iterator.__next__)
        util_logger.debug("Fetching next iterator element, eagerly")

    def __next__(self) -> _T:
        try:
            assert self._future
            res = self._future.result()
            self._begin_next()
            return res
        except StopIteration:
            self._cleanup()
            raise

    def _cleanup(self) -> None:
        util_logger.debug("Cleaning up eager iterator")
        if self._own_pool:
            self._pool.shutdown()

    def __del__(self) -> None:
        # Ensure the threadpool is cleaned up in the case where the
        # iterator is not exhausted. For more information on __del__:
        # https://docs.python.org/3/reference/datamodel.html#object.__del__
        self._cleanup()
        super_del = getattr(super(), "__del__", lambda: None)
        super_del()


class _EagerBufferedIterator(Iterator[_T]):
    def __init__(
        self,
        iterator: Iterator[_T],
        max_pending: int = 1,
        pool: futures.Executor | None = None,
    ):
        super().__init__()
        self.iterator = iterator
        self.max_pending = max_pending
        self._pool = pool or futures.ThreadPoolExecutor()
        self._own_pool = pool is None
        self._pending_results: deque[futures.Future[_T]] = deque()
        self._lock = threading.Lock()
        self._begin_next()

    def __next__(self) -> _T:
        try:
            res = self._pending_results[0].result()
            self._pending_results.popleft()
            self._begin_next()
            return res
        except StopIteration:
            self._cleanup()
            raise

    def _begin_next(self) -> None:
        def _fut_done(fut: futures.Future[_T]) -> None:
            util_logger.debug("Finished fetching next iterator element, eagerly")
            if fut.exception() is None:
                self._begin_next()

        with self._lock:
            not_running = len(self._pending_results) == 0 or self._pending_results[-1].done()
            if len(self._pending_results) < self.max_pending and not_running:
                _future = self._pool.submit(self.iterator.__next__)
                util_logger.debug("Fetching next iterator element, eagerly")
                _future.add_done_callback(_fut_done)
                self._pending_results.append(_future)
            assert len(self._pending_results) <= self.max_pending

    def _cleanup(self) -> None:
        util_logger.debug("Cleaning up eager iterator")
        if self._own_pool:
            self._pool.shutdown()

    def __del__(self) -> None:
        # Ensure the threadpool is cleaned up in the case where the
        # iterator is not exhausted. For more information on __del__:
        # https://docs.python.org/3/reference/datamodel.html#object.__del__
        self._cleanup()
        super_del = getattr(super(), "__del__", lambda: None)
        super_del()


class Encoder(abc.ABC):
    """Base class for obs encoders.

    To define a custom encoder, two methods must be implemented:

    - ``register``: defines how the encoder will be fitted to the data.
    - ``transform``: defines how the encoder will be applied to the data
      in order to create an obs_tensor.

    See the implementation of ``DefaultEncoder`` for an example.
    """

    @abc.abstractmethod
    def register(self, obs: pd.DataFrame) -> None:
        """Register the encoder with obs."""
        pass

    @abc.abstractmethod
    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform the obs DataFrame into a DataFrame of encoded values."""
        pass

    @property
    def name(self) -> str:
        return self.__class__.__name__


class DefaultEncoder(Encoder):
    """Default encoder based on LabelEncoder."""

    def __init__(self, col: str) -> None:
        self._encoder = LabelEncoder()
        self.col = col

    def register(self, obs: pd.DataFrame) -> None:
        self._encoder.fit(obs[self.col].unique())

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        return self._encoder.transform(df[self.col])  # type: ignore

    @property
    def name(self) -> str:
        return self.col

    @property
    def classes_(self):  # type: ignore
        return self._encoder.classes_


@define
class _SOMAChunk:
    """Return type of ``_ObsAndXSOMAIterator`` that pairs a chunk of ``obs`` rows with the

    respective rows from the ``X`` matrix.

    Lifecycle:
        experimental
    """

    obs: pd.DataFrame
    X: scipy.sparse.spmatrix
    stats: Stats

    def __len__(self) -> int:
        return len(self.obs)


Encoders = dict[str, LabelEncoder]
"""A dictionary of ``LabelEncoder``s keyed by the ``obs`` column name."""


@define
class Stats:
    """Statistics about the data retrieved by ``ExperimentDataPipe`` via SOMA API. This is useful

    for assessing the read throughput of SOMA data.

    Lifecycle:
        experimental
    """

    n_obs: int = 0
    """The total number of obs rows retrieved"""

    nnz: int = 0
    """The total number of values retrieved"""

    elapsed: int = 0
    """The total elapsed time in seconds for retrieving all batches"""

    n_soma_chunks: int = 0
    """The number of chunks retrieved"""

    def __str__(self) -> str:
        return (
            f"{self.n_soma_chunks=}, {self.n_obs=}, {self.nnz=}, "
            f"elapsed={timedelta(seconds=self.elapsed)}"
        )

    def __add__(self, other: Stats) -> Stats:
        self.n_obs += other.n_obs
        self.nnz += other.nnz
        self.elapsed += other.elapsed
        self.n_soma_chunks += other.n_soma_chunks
        return self


@contextmanager
def _open_experiment(
    uri: str,
    aws_region: str | None = None,
) -> soma.Experiment:
    """Internal method for opening a SOMA ``Experiment`` as a context manager."""
    context = get_default_soma_context().replace(
        tiledb_config={"vfs.s3.region": aws_region} if aws_region else {}
    )

    with soma.Experiment.open(uri, context=context) as exp:
        yield exp


class _ObsAndXSOMAIterator(Iterator[_SOMAChunk]):
    """Iterates the SOMA chunks of corresponding ``obs`` and ``X`` data. This is an internal class,

    not intended for public use.
    """

    X: soma.SparseNDArray
    """A handle to the full X data of the SOMA ``Experiment``"""

    obs_joinids_chunks_iter: Iterator[npt.NDArray[np.int64]]

    var_joinids: npt.NDArray[np.int64]
    """The ``var`` joinids to be retrieved from the SOMA ``Experiment``"""

    def __init__(
        self,
        obs: soma.DataFrame,
        X: soma.SparseNDArray,
        obs_column_names: Sequence[str],
        obs_joinids_chunked: list[npt.NDArray[np.int64]],
        var_joinids: npt.NDArray[np.int64],
        shuffle_chunk_count: int | None = None,
        shuffle_rng: Generator | None = None,
    ):
        self.obs = obs
        self.X = X
        self.obs_column_names = obs_column_names
        if shuffle_chunk_count:
            assert shuffle_rng is not None

            # At the start of this step, `obs_joinids_chunked` is a list of one dimensional
            # numpy arrays. Each numpy array corresponds to a chunk of contiguous rows in `obs`.
            # Critically, `obs_joinids_chunked` is randomly ordered where each chunk is
            # from a random section of `obs`.
            # We then take `shuffle_chunk_count` of these in order, concatenate them into
            # a larger numpy array and shuffle this larger numpy array.
            # The result is again a list of numpy arrays.
            self.obs_joinids_chunks_iter = (
                shuffle_rng.permutation(np.concatenate(grouped_chunks))
                for grouped_chunks in list_split(obs_joinids_chunked, shuffle_chunk_count)
            )
        else:
            self.obs_joinids_chunks_iter = iter(obs_joinids_chunked)
        self.var_joinids = var_joinids
        self.shuffle_chunk_count = shuffle_chunk_count

    def __next__(self) -> _SOMAChunk:
        pytorch_logger.debug("Retrieving next SOMA chunk...")
        start_time = time()

        # If no more chunks to iterate through, raise StopIteration, as all iterators
        # do when at end
        obs_joinids_chunk = next(self.obs_joinids_chunks_iter)

        obs_batch = (
            self.obs.read(
                coords=(obs_joinids_chunk,),
                column_names=self.obs_column_names,
            )
            .concat()
            .to_pandas()
            .set_index("soma_joinid")
        )
        assert obs_batch.shape[0] == obs_joinids_chunk.shape[0]

        # handle case of empty result (first batch has 0 rows)
        if len(obs_batch) == 0:
            raise StopIteration

        # reorder obs rows to match obs_joinids_chunk ordering, which may be shuffled
        obs_batch = obs_batch.reindex(obs_joinids_chunk, copy=False)

        # note: the `blockwise` call is employed for its ability to reindex the axes of the sparse
        # matrix, but the blockwise iteration feature is not used (block_size is set to retrieve
        # the chunk as a single block)
        scipy_iter = (
            self.X.read(coords=(obs_joinids_chunk, self.var_joinids))
            .blockwise(axis=0, size=len(obs_joinids_chunk), eager=False)
            .scipy(compress=True)
        )
        X_batch, _ = next(scipy_iter)
        assert obs_batch.shape[0] == X_batch.shape[0]

        stats = Stats()
        stats.n_obs += X_batch.shape[0]
        stats.nnz += X_batch.nnz
        stats.elapsed += int(time() - start_time)
        stats.n_soma_chunks += 1

        pytorch_logger.debug(f"Retrieved SOMA chunk: {stats}")
        return _SOMAChunk(obs=obs_batch, X=X_batch, stats=stats)


def list_split(arr_list: list[Any], sublist_len: int) -> list[list[Any]]:
    """Splits a python list into a list of sublists where each sublist is of size `sublist_len`.

    TODO: Replace with `itertools.batched` when Python 3.12 becomes the minimum supported version.
    """
    i = 0
    result = []
    while i < len(arr_list):
        if (i + sublist_len) >= len(arr_list):
            result.append(arr_list[i:])
        else:
            result.append(arr_list[i : i + sublist_len])

        i += sublist_len

    return result


def run_gc() -> tuple[tuple[Any, Any, Any], tuple[Any, Any, Any]]:
    proc = psutil.Process(os.getpid())

    pre_gc = proc.memory_full_info(), psutil.virtual_memory(), psutil.swap_memory()
    gc.collect()
    post_gc = proc.memory_full_info(), psutil.virtual_memory(), psutil.swap_memory()

    pytorch_logger.debug(f"gc:  pre={pre_gc}")
    pytorch_logger.debug(f"gc: post={post_gc}")

    return pre_gc, post_gc


class _ObsAndXIterator(Iterator[ObsAndXDatum]):
    """Iterates through a set of ``obs`` and corresponding ``X`` rows, where the rows to be

    returned are specified by the ``obs_tables_iter`` argument. For the specified ``obs` rows,
    the corresponding ``X`` data is loaded and joined together. It is returned from this iterator
    as 2-tuples of ``X`` and obs Tensors.

    Internally manages the retrieval of data in SOMA-sized chunks, fetching the next chunk of SOMA
    data as needed. Supports fetching the data in an eager manner, where the next SOMA chunk is
    fetched while the current chunk is being read. This is an internal class, not intended for
    public use.
    """

    soma_chunk_iter: _SOMAChunk | None
    """The iterator for SOMA chunks of paired obs and X data"""

    soma_chunk: _SOMAChunk | None
    """The current SOMA chunk of obs and X data"""

    i: int = -1
    """Index into current obs ``SOMA`` chunk"""

    def __init__(
        self,
        obs: soma.DataFrame,
        X: soma.SparseNDArray,
        obs_column_names: Sequence[str],
        obs_joinids_chunked: list[npt.NDArray[np.int64]],
        var_joinids: npt.NDArray[np.int64],
        batch_size: int,
        encoders: list[Encoder],
        stats: Stats,
        return_sparse_X: bool,
        use_eager_fetch: bool,
        shuffle_chunk_count: int | None = None,
        shuffle_rng: Generator | None = None,
    ) -> None:
        self.soma_chunk_iter = _ObsAndXSOMAIterator(
            obs,
            X,
            obs_column_names,
            obs_joinids_chunked,
            var_joinids,
            shuffle_chunk_count,
            shuffle_rng,
        )
        if use_eager_fetch:
            self.soma_chunk_iter = _EagerIterator(self.soma_chunk_iter)
        self.soma_chunk = None
        self.var_joinids = var_joinids
        self.batch_size = batch_size
        self.return_sparse_X = return_sparse_X
        self.encoders = encoders
        self.stats = stats
        self.max_process_mem_usage_bytes = 0
        self.X_dtype = X.schema[2].type.to_pandas_dtype()

    def __next__(self) -> ObsAndXDatum:
        """Read the next torch batch, possibly across multiple soma chunks."""
        obs: pd.DataFrame = pd.DataFrame()
        X: sparse.csr_matrix = sparse.csr_matrix((0, len(self.var_joinids)), dtype=self.X_dtype)

        while len(obs) < self.batch_size:
            try:
                obs_partial, X_partial = self._read_partial_torch_batch(self.batch_size - len(obs))
                obs = pd.concat([obs, obs_partial], axis=0)
                X = sparse.vstack([X, X_partial])
            except StopIteration:
                break

        if len(obs) == 0:
            raise StopIteration

        obs_encoded = pd.DataFrame()

        for enc in self.encoders:
            obs_encoded[enc.name] = enc.transform(obs)

        # `to_numpy()` avoids copying the numpy array data
        obs_tensor = torch.from_numpy(obs_encoded.to_numpy())

        if not self.return_sparse_X:
            X_tensor = torch.from_numpy(X.todense())
        else:
            coo = X.tocoo()

            X_tensor = torch.sparse_coo_tensor(
                # Note: The `np.array` seems unnecessary, but PyTorch warns bare array
                # is "extremely slow"
                indices=torch.from_numpy(np.array([coo.row, coo.col])),
                values=coo.data,
                size=coo.shape,
            )

        if self.batch_size == 1:
            X_tensor = X_tensor[0]
            obs_tensor = obs_tensor[0]

        return X_tensor, obs_tensor

    def _read_partial_torch_batch(self, batch_size: int) -> ObsAndXDatum:
        """Reads a torch-size batch of data from the current SOMA chunk, returning a torch-size

        batch whose size may contain fewer rows than the requested ``batch_size``. This can happen
        when the remaining rows in the current SOMA chunk are fewer than the requested
        ``batch_size``.
        """
        if self.soma_chunk is None or not (0 <= self.i < len(self.soma_chunk)):
            # GC memory from previous soma_chunk
            self.soma_chunk = None
            mem_info = run_gc()
            self.max_process_mem_usage_bytes = max(
                self.max_process_mem_usage_bytes, mem_info[0][0].uss
            )

            self.soma_chunk: _SOMAChunk = next(self.soma_chunk_iter)
            self.stats += self.soma_chunk.stats
            self.i = 0

            pytorch_logger.debug(f"Retrieved SOMA chunk totals: {self.stats}")

        obs_batch = self.soma_chunk.obs
        X_batch = self.soma_chunk.X

        safe_batch_size = min(batch_size, len(obs_batch) - self.i)
        slice_ = slice(self.i, self.i + safe_batch_size)
        assert slice_.stop <= obs_batch.shape[0]

        obs_rows = obs_batch.iloc[slice_]
        assert obs_rows.index.is_unique
        assert safe_batch_size == obs_rows.shape[0]

        X_csr_scipy = X_batch[slice_]
        assert obs_rows.shape[0] == X_csr_scipy.shape[0]

        self.i += safe_batch_size

        return obs_rows, X_csr_scipy


class ExperimentDataPipe(pipes.IterDataPipe[Dataset[ObsAndXDatum]]):  # type: ignore
    r"""An :class:`torchdata.datapipes.iter.IterDataPipe` that reads ``obs`` and ``X`` data from a

    :class:`tiledbsoma.Experiment`, based upon the specified queries along the ``obs`` and ``var``
    axes. Provides an iterator over these data when the object is passed to Python's built-in
    ``iter`` function.

    >>> for batch in iter(ExperimentDataPipe(...)):
            X_batch, y_batch = batch

    The ``batch_size`` parameter controls the number of rows of ``obs`` and ``X`` data that are
    returned in each iteration. If the ``batch_size`` is 1, then each Tensor will have rank 1:

    >>> (tensor([0., 0., 0., 0., 0., 1., 0., 0., 0.]),  # X data
         tensor([2415,    0,    0], dtype=torch.int64)) # obs data, encoded

    For larger ``batch_size`` values, the returned Tensors will have rank 2:

    >>> DataLoader(..., batch_size=3, ...):
        (tensor([[0., 0., 0., 0., 0., 1., 0., 0., 0.],     # X batch
                 [0., 0., 0., 0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0., 0., 0., 0.]]),
         tensor([[2415,    0,    0],                       # obs batch
                 [2416,    0,    4],
                 [2417,    0,    3]], dtype=torch.int64))

    The ``return_sparse_X`` parameter controls whether the ``X`` data is returned as a dense or
    sparse :class:`torch.Tensor`. If the model supports use of sparse :class:`torch.Tensor`\ s,
    this will reduce memory usage.

    The ``obs_column_names`` parameter determines the data columns that are returned in the
    ``obs`` Tensor. The first element is always the ``soma_joinid`` of the ``obs``
    :class:`pandas.DataFrame` (or, equivalently, the ``soma_dim_0`` of the ``X`` matrix).
    The remaining elements are the ``obs`` columns specified by ``obs_column_names``,
    and string-typed columns are encoded as integer values. If needed, these values can be decoded
    by obtaining the encoder for a given ``obs`` column name and calling its ``inverse_transform``
    method:

    >>> exp_data_pipe.obs_encoders["<obs_attr_name>"].inverse_transform(encoded_values)

    Lifecycle:
        experimental
    """

    _initialized: bool

    _obs_joinids: npt.NDArray[np.int64] | None

    _var_joinids: npt.NDArray[np.int64] | None

    _encoders: list[Encoder]

    _stats: Stats

    _shuffle_rng: Generator | None

    # TODO: Consider adding another convenience method wrapper to construct this object whose
    #  signature is more closely aligned with get_anndata() params
    #  (i.e. "exploded" AxisQuery params).
    def __init__(
        self,
        experiment: soma.Experiment,
        measurement_name: str = "RNA",
        X_name: str = "raw",
        obs_query: soma.AxisQuery | None = None,
        var_query: soma.AxisQuery | None = None,
        obs_column_names: Sequence[str] = (),
        batch_size: int = 1,
        shuffle: bool = True,
        seed: int | None = None,
        return_sparse_X: bool = False,
        soma_chunk_size: int | None = 64,
        use_eager_fetch: bool = True,
        encoders: list[Encoder] | None = None,
        shuffle_chunk_count: int | None = 2000,
    ) -> None:
        r"""Construct a new ``ExperimentDataPipe``.

        Args:
            experiment:
                The :class:`tiledbsoma.Experiment` from which to read data.
            measurement_name:
                The name of the :class:`tiledbsoma.Measurement` to read. Defaults to ``"RNA"``.
            X_name:
                The name of the X layer to read. Defaults to ``"raw"``.
            obs_query:
                The query used to filter along the ``obs`` axis. If not specified, all ``obs`` and
                ``X`` data will be returned, which can be very large.
            var_query:
                The query used to filter along the ``var`` axis. If not specified, all ``var``
                columns (genes/features) will be returned.
            obs_column_names:
                The names of the ``obs`` columns to return. The ``soma_joinid`` index "column" does
                not need to be specified and will always be returned. If not specified, only the
                ``soma_joinid`` will be returned.
            batch_size:
                The number of rows of ``obs`` and ``X`` data to return in each iteration. Defaults
                to ``1``. A value of ``1`` will result in :class:`torch.Tensor` of rank 1 being
                returns (a single row); larger values will result in :class:`torch.Tensor`\ s of
                rank 2 (multiple rows).
            shuffle:
                Whether to shuffle the ``obs`` and ``X`` data being returned. Defaults to ``True``.
                For performance reasons, shuffling is not performed globally across all rows, but
                rather in chunks. More specifically, we select ``shuffle_chunk_count``
                non-contiguous chunks across all the observations
                in the query, concatenate the chunks and shuffle the associated observations.
                The randomness of the shuffling is therefore determined by the
                (``soma_chunk_size``, ``shuffle_chunk_count``) selection. The default values have
                been determined to yield a good trade-off between randomness and performance.
                Further tuning may be required for different type of models. Note that memory usage
                is correlated to the product ``soma_chunk_size * shuffle_chunk_count``.
            seed:
                The random seed used for shuffling. Defaults to ``None`` (no seed). This *must* be
                specified when using :class:`torch.nn.parallel.DistributedDataParallel` to ensure
                data partitions are disjoint across worker processes.
            return_sparse_X:
                Controls whether the ``X`` data is returned as a dense or sparse
                :class:`torch.Tensor`. As ``X`` data is very sparse, setting this to ``True`` will
                reduce memory usage, if the model supports use of sparse :class:`torch.Tensor`\ s.
                Defaults to ``False``, since sparse :class:`torch.Tensor`\ s are still experimental
                in PyTorch.
            soma_chunk_size:
                The number of ``obs``/``X`` rows to retrieve when reading data from SOMA. This
                impacts two aspects of this class's behavior: 1) The maximum memory utilization,
                with larger values providing better read performance, but also requiring more
                memory; 2) The granularity of the global shuffling step (see ``shuffle`` parameter
                for details). The default value of 64 works well in conjunction with the default
                ``shuffle_chunk_count`` value.
            use_eager_fetch:
                Fetch the next SOMA chunk of ``obs`` and ``X`` data immediately after a previously
                fetched SOMA chunk is made available for processing via the iterator. This allows
                network (or filesystem) requests to be made in parallel with client-side processing
                 of the SOMA data, potentially improving overall performance at the cost of
                 doubling memory utilization. Defaults to ``True``.
            shuffle_chunk_count:
                The number of contiguous blocks (chunks) of rows sampled to then concatenate
                and shuffle. Larger numbers correspond to more randomness per training batch.
                If ``shuffle == False``, this parameter is ignored. Defaults to ``2000``.
            encoders:
                Specify custom encoders to be used. If not specified, a LabelEncoder will be
                created and used for each column in ``obs_column_names``. If specified, only
                columns for which an encoder has been registered will be returned in the
                ``obs`` tensor.

        Lifecycle:
            experimental
        """
        self.exp_uri = experiment.uri
        self.aws_region = experiment.context.tiledb_ctx.config().get("vfs.s3.region")
        self.measurement_name = measurement_name
        self.layer_name = X_name
        self.obs_query = obs_query
        self.var_query = var_query
        self.obs_column_names = obs_column_names
        self.batch_size = batch_size
        self.return_sparse_X = return_sparse_X
        self.soma_chunk_size = soma_chunk_size
        self.use_eager_fetch = use_eager_fetch
        self._stats = Stats()
        self._custom_encoders = encoders
        self._encoders = []
        self._obs_joinids = None
        self._var_joinids = None
        self._shuffle_chunk_count = shuffle_chunk_count if shuffle else None
        self._shuffle_rng = np.random.default_rng(seed) if shuffle else None
        self._initialized = False

        if "soma_joinid" not in self.obs_column_names:
            self.obs_column_names = ["soma_joinid", *self.obs_column_names]

    def _init(self) -> None:
        if self._initialized:
            return

        pytorch_logger.debug("Initializing ExperimentDataPipe")

        with _open_experiment(self.exp_uri, self.aws_region) as exp:
            query = exp.axis_query(
                measurement_name=self.measurement_name,
                obs_query=self.obs_query,
                var_query=self.var_query,
            )

            # The to_numpy() call is a workaround for a possible bug in TileDB-SOMA:
            # https://github.com/single-cell-data/TileDB-SOMA/issues/1456
            self._obs_joinids = query.obs_joinids().to_numpy()
            self._var_joinids = query.var_joinids().to_numpy()

            self._encoders = self._build_obs_encoders(query)

        self._initialized = True

    @staticmethod
    def _subset_ids_to_partition(
        ids_chunked: list[npt.NDArray[np.int64]],
        partition_index: int,
        num_partitions: int,
    ) -> list[npt.NDArray[np.int64]]:
        """Returns a single partition of the obs_joinids_chunked (a 2D ndarray),

        based upon the current process's distributed rank and world size.
        """
        # subset to a single partition
        # typing does not reflect that is actually a list of 2D NDArrays
        partition_indices = np.array_split(range(len(ids_chunked)), num_partitions)
        partition = [ids_chunked[i] for i in partition_indices[partition_index]]

        if pytorch_logger.isEnabledFor(logging.DEBUG) and len(partition) > 0:
            pytorch_logger.debug(
                f"Process {os.getpid()} handling partition {partition_index + 1} "
                f"of {num_partitions}, partition_size={sum([len(chunk) for chunk in partition])}"
            )

        return partition

    @staticmethod
    def _compute_partitions(
        loader_partition: int,
        loader_partitions: int,
        dist_partition: int,
        num_dist_partitions: int,
    ) -> tuple[int, int]:
        # NOTE: Can alternately use a `worker_init_fn` to split among workers split workload
        total_partitions = num_dist_partitions * loader_partitions
        partition = dist_partition * loader_partitions + loader_partition
        return partition, total_partitions

    def __iter__(self) -> Iterator[ObsAndXDatum]:
        self._init()
        assert self._obs_joinids is not None
        assert self._var_joinids is not None

        if self.soma_chunk_size is None:
            # set soma_chunk_size to utilize ~1 GiB of RAM per SOMA chunk; assumes 95% X data
            # sparsity, 8 bytes for the X value and 8 bytes for the sparse matrix indices,
            # and a 100% working memory overhead (2x).
            X_row_memory_size = 0.05 * len(self._var_joinids) * 8 * 3 * 2
            self.soma_chunk_size = int((1 * 1024**3) / X_row_memory_size)
        pytorch_logger.debug(f"Using {self.soma_chunk_size=}")

        if (
            self.return_sparse_X
            and torch.utils.data.get_worker_info()
            and torch.utils.data.get_worker_info().num_workers > 0
        ):
            raise NotImplementedError(
                "torch does not work with sparse tensors in multi-processing mode "
                "(see https://github.com/pytorch/pytorch/issues/20248)"
            )

        # chunk the obs joinids into batches of size soma_chunk_size
        obs_joinids_chunked = self._chunk_ids(self._obs_joinids, self.soma_chunk_size)

        # globally shuffle the chunks, if requested
        if self._shuffle_rng:
            self._shuffle_rng.shuffle(obs_joinids_chunked)

        # subset to a single partition, as needed for distributed training and multi-processing
        # data loading
        worker_info = torch.utils.data.get_worker_info()
        partition, partitions = self._compute_partitions(
            loader_partition=worker_info.id if worker_info else 0,
            loader_partitions=worker_info.num_workers if worker_info else 1,
            dist_partition=dist.get_rank() if dist.is_initialized() else 0,
            num_dist_partitions=dist.get_world_size() if dist.is_initialized() else 1,
        )
        obs_joinids_chunked_partition: list[npt.NDArray[np.int64]] = self._subset_ids_to_partition(
            obs_joinids_chunked, partition, partitions
        )

        with _open_experiment(self.exp_uri, self.aws_region) as exp:
            obs_and_x_iter = _ObsAndXIterator(
                obs=exp.obs,
                X=exp.ms[self.measurement_name].X[self.layer_name],
                obs_column_names=self.obs_column_names,
                obs_joinids_chunked=obs_joinids_chunked_partition,
                var_joinids=self._var_joinids,
                batch_size=self.batch_size,
                encoders=self._encoders,
                stats=self._stats,
                return_sparse_X=self.return_sparse_X,
                use_eager_fetch=self.use_eager_fetch,
                shuffle_rng=self._shuffle_rng,
                shuffle_chunk_count=self._shuffle_chunk_count,
            )

            yield from obs_and_x_iter

            pytorch_logger.debug(
                "max process memory usage="
                f"{obs_and_x_iter.max_process_mem_usage_bytes / (1024 ** 3):.3f} GiB"
            )

    @staticmethod
    def _chunk_ids(ids: npt.NDArray[np.int64], chunk_size: int) -> list[npt.NDArray[np.int64]]:
        num_chunks = max(1, ceil(len(ids) / chunk_size))
        pytorch_logger.debug(
            f"Shuffling {len(ids)} obs joinids into {num_chunks} chunks of {chunk_size}"
        )
        return np.array_split(ids, num_chunks)

    def __len__(self) -> int:
        self._init()
        assert self._obs_joinids is not None

        return len(self._obs_joinids)

    def __getitem__(self, index: int) -> ObsAndXDatum:
        raise NotImplementedError("IterDataPipe can only be iterated")

    def _build_obs_encoders(self, query: soma.ExperimentAxisQuery) -> list[Encoder]:
        pytorch_logger.debug("Initializing encoders")

        encoders = []
        obs = query.obs(column_names=self.obs_column_names).concat().to_pandas()

        if self._custom_encoders:
            # Register all the custom encoders with obs
            for enc in self._custom_encoders:
                enc.register(obs)
                encoders.append(enc)
        else:
            # Create one DefaultEncoder for each column, and register it with obs
            for col in self.obs_column_names:
                if obs[col].dtype in [object]:
                    enc = DefaultEncoder(col)
                    enc.register(obs)
                    encoders.append(enc)

        return encoders

    # TODO: This does not work in multiprocessing mode, as child process's stats are not collected
    def stats(self) -> Stats:
        """Get data loading stats for this

        :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`.

        Returns
        -------
            The :class:`cellxgene_census.experimental.ml.pytorch.Stats` object for this
            :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`.

        Lifecycle:
            experimental
        """
        return self._stats

    @property
    def shape(self) -> tuple[int, int]:
        """Get the shape of the data that will be returned by this

        :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`.
        This is the number of obs (cell) and var (feature) counts in the returned data. If used in
        multiprocessing mode (i.e. :class:`torch.utils.data.DataLoader`
        instantiated with num_workers > 0), the obs (cell) count will reflect
        the size of the partition of the data assigned to the active process.

        Returns
        -------
            A 2-tuple of ``int``s, for obs and var counts, respectively.

        Lifecycle:
            experimental
        """
        self._init()
        assert self._obs_joinids is not None
        assert self._var_joinids is not None

        return len(self._obs_joinids), len(self._var_joinids)

    @property
    def obs_encoders(self) -> Encoders:
        """Returns a dictionary of :class:`sklearn.preprocessing.LabelEncoder` objects, keyed on

        ``obs`` column names, which were used to encode the ``obs`` column values.

        These encoders can be used to decode the encoded values as follows:

        >>> exp_data_pipe.obs_encoders["<obs_attr_name>"].inverse_transform(encoded_values)

        Returns
        -------
            A ``dict[str, LabelEncoder]``, mapping column names to :class:`sklearn.preprocessing.
            LabelEncoder` objects.
        """
        self._init()
        assert self._encoders is not None

        return {enc.name: enc for enc in self._encoders}


# Note: must be a top-level function (and not a lambda), to play nice with multiprocessing pickling
def _collate_noop(x: Any) -> Any:
    return x


# TODO: Move into somacore.ExperimentAxisQuery
def experiment_dataloader(
    datapipe: pipes.IterDataPipe,
    num_workers: int = 0,
    **dataloader_kwargs: Any,
) -> DataLoader:
    """Factory method for :class:`torch.utils.data.DataLoader`. This method can be used to safely

    instantiate a :class:`torch.utils.data.DataLoader` that works with
    :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`, since some of the
    :class:`torch.utils.data.DataLoader` constructor parameters are not applicable when using a
    :class:`torchdata.datapipes.iter.IterDataPipe` (``shuffle``, ``batch_size``, ``sampler``,
    ``batch_sampler``,``collate_fn``).

    Args:
        datapipe:
            An :class:`torchdata.datapipes.iter.IterDataPipe`, which can be an
            :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe` or any other
            :class:`torchdata.datapipes.iter.IterDataPipe` that has been chained to the
            :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`.
        num_workers:
            Number of worker processes to use for data loading. If ``0``, data will be loaded in
            the main process.
        **dataloader_kwargs:
            Additional keyword arguments to pass to the :class:`torch.utils.data.DataLoader`
            constructor, except for ``shuffle``, ``batch_size``, ``sampler``, ``batch_sampler``,
            and ``collate_fn``, which are not supported when using
            :class:`cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`.

    Returns
    -------
        A :class:`torch.utils.data.DataLoader`.

    Raises
    ------
        ValueError: if any of the ``shuffle``, ``batch_size``, ``sampler``, ``batch_sampler``,
        or ``collate_fn`` params are passed as keyword arguments.

    Lifecycle:
        experimental
    """
    unsupported_dataloader_args = [
        "shuffle",
        "batch_size",
        "sampler",
        "batch_sampler",
        "collate_fn",
    ]
    if set(unsupported_dataloader_args).intersection(dataloader_kwargs.keys()):
        raise ValueError(
            f"The {','.join(unsupported_dataloader_args)} DataLoader params are not supported"
        )

    if num_workers > 0:
        _init_multiprocessing()

    return DataLoader(
        datapipe,
        batch_size=None,  # batching is handled by our ExperimentDataPipe
        num_workers=num_workers,
        # avoid use of default collator, which adds an extra (3rd) dimension to the tensor batches
        collate_fn=_collate_noop,
        # shuffling is handled by our ExperimentDataPipe
        shuffle=False,
        **dataloader_kwargs,
    )


def _init_multiprocessing() -> None:
    """Ensures use of "spawn" for starting child processes with multiprocessing.

    Forked processes are known to be problematic:
      https://pytorch.org/docs/stable/notes/multiprocessing.html#avoiding-and-fighting-deadlocks
    Also, CUDA does not support forked child processes:
      https://pytorch.org/docs/stable/notes/multiprocessing.html#cuda-in-multiprocessing

    """
    torch.multiprocessing.set_start_method("fork", force=True)
    orig_start_method = torch.multiprocessing.get_start_method()
    if orig_start_method != "spawn":
        if orig_start_method:
            pytorch_logger.warning(
                "switching torch multiprocessing start method from "
                f'"{torch.multiprocessing.get_start_method()}" to "spawn"'
            )
        torch.multiprocessing.set_start_method("spawn", force=True)


class BatchEncoder(Encoder):
    """Concatenates and encodes several columns."""

    def __init__(self, cols: list[str]):
        self.cols = cols
        from sklearn.preprocessing import LabelEncoder

        self._encoder = LabelEncoder()

    def transform(self, df: pd.DataFrame):
        import functools

        arr = functools.reduce(lambda a, b: a + b, [df[c].astype(str) for c in self.cols])
        return self._encoder.transform(arr)

    def register(self, obs: pd.DataFrame):
        import functools

        arr = functools.reduce(lambda a, b: a + b, [obs[c].astype(str) for c in self.cols])
        self._encoder.fit(arr.unique())

    @property
    def name(self) -> str:
        return "batch"

    @property
    def classes_(self):
        return self._encoder.classes_


class CensusSCVIDataModule(LightningDataModule):
    """Lightning data module for CxG Census.

    Parameters
    ----------
    *args
        Positional arguments passed to
        :class:`~cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`.
    batch_keys
        List of obs column names concatenated to form the batch column.
    train_size
        Fraction of data to use for training.
    split_seed
        Seed for data split.
    dataloader_kwargs
        Keyword arguments passed into
        :func:`~cellxgene_census.experimental.ml.pytorch.experiment_dataloader`.
    **kwargs
        Additional keyword arguments passed into
        :class:`~cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`. Must not include
        ``obs_column_names``.
    """

    _TRAIN_KEY = "train"
    _VALIDATION_KEY = "validation"

    def __init__(
        self,
        *args,
        batch_keys: list[str] | None = None,
        train_size: float | None = None,
        split_seed: int | None = None,
        dataloader_kwargs: dict[str, any] | None = None,
        **kwargs,
    ):
        super().__init__()
        self.datapipe_args = args
        self.datapipe_kwargs = kwargs
        self.batch_keys = batch_keys
        self.train_size = train_size
        self.split_seed = split_seed
        self.dataloader_kwargs = dataloader_kwargs or {}

    @property
    def batch_keys(self) -> list[str]:
        """List of obs column names concatenated to form the batch column."""
        if not hasattr(self, "_batch_keys"):
            raise AttributeError("`batch_keys` not set.")
        return self._batch_keys

    @batch_keys.setter
    def batch_keys(self, value: list[str] | None):
        if value is None or not isinstance(value, list):
            raise ValueError("`batch_keys` must be a list of strings.")
        self._batch_keys = value

    @property
    def obs_column_names(self) -> list[str]:
        """Passed to :class:`~cellxgene_census.experimental.ml.pytorch.ExperimentDataPipe`."""
        if hasattr(self, "_obs_column_names"):
            return self._obs_column_names

        obs_column_names = []
        if self.batch_keys is not None:
            obs_column_names.extend(self.batch_keys)

        self._obs_column_names = obs_column_names
        return self._obs_column_names

    @property
    def split_seed(self) -> int:
        """Seed for data split."""
        if not hasattr(self, "_split_seed"):
            raise AttributeError("`split_seed` not set.")
        return self._split_seed

    @split_seed.setter
    def split_seed(self, value: int | None):
        if value is not None and not isinstance(value, int):
            raise ValueError("`split_seed` must be an integer.")
        self._split_seed = value or 0

    @property
    def train_size(self) -> float:
        """Fraction of data to use for training."""
        if not hasattr(self, "_train_size"):
            raise AttributeError("`train_size` not set.")
        return self._train_size

    @train_size.setter
    def train_size(self, value: float | None):
        if value is not None and not isinstance(value, float):
            raise ValueError("`train_size` must be a float.")
        elif value is not None and (value < 0.0 or value > 1.0):
            raise ValueError("`train_size` must be between 0.0 and 1.0.")
        self._train_size = value or 1.0

    @property
    def validation_size(self) -> float:
        """Fraction of data to use for validation."""
        if not hasattr(self, "_train_size"):
            raise AttributeError("`validation_size` not available.")
        return 1.0 - self.train_size

    @property
    def weights(self) -> dict[str, float]:
        """Passed to :meth:`~cellxgene_census.experimental.ml.ExperimentDataPipe.random_split`."""
        if not hasattr(self, "_weights"):
            self._weights = {self._TRAIN_KEY: self.train_size}
            if self.validation_size > 0.0:
                self._weights[self._VALIDATION_KEY] = self.validation_size
        return self._weights

    @property
    def datapipe(self) -> ExperimentDataPipe:
        """Experiment data pipe."""
        if not hasattr(self, "_datapipe"):
            encoder = BatchEncoder(self.obs_column_names)
            self._datapipe = ExperimentDataPipe(
                *self.datapipe_args,
                obs_column_names=self.obs_column_names,
                encoders=[encoder],
                **self.datapipe_kwargs,
            )
        return self._datapipe

    def setup(self, stage: str | None = None):
        """Set up the train and validation data pipes."""
        datapipes = self.datapipe.random_split(weights=self.weights, seed=self.split_seed)
        self._train_datapipe = datapipes[0]
        if self.validation_size > 0.0:
            self._validation_datapipe = datapipes[1]
        else:
            self._validation_datapipe = None

    def train_dataloader(self):
        """Training data loader."""
        return experiment_dataloader(self._train_datapipe, **self.dataloader_kwargs)

    def val_dataloader(self):
        """Validation data loader."""
        if self._validation_datapipe is not None:
            return experiment_dataloader(self._validation_datapipe, **self.dataloader_kwargs)

    @property
    def n_obs(self) -> int:
        """Number of observations in the query.

        Necessary in scvi-tools to compute a heuristic of ``max_epochs``.
        """
        return self.datapipe.shape[0]

    @property
    def n_vars(self) -> int:
        """Number of features in the query.

        Necessary in scvi-tools to initialize the actual layers in the model.

        """
        return self.datapipe.shape[1]

    @property
    def n_batch(self) -> int:
        """
        Number of unique batches (after concatenation of ``batch_keys``). Necessary in scvi-tools

        so that the model knows how to one-hot encode batches.

        """
        return self.get_n_classes("batch")

    def get_n_classes(self, key: str) -> int:
        """Return the number of classes for a given obs column."""
        return len(self.datapipe.obs_encoders[key].classes_)

    def on_before_batch_transfer(
        self,
        batch: tuple[torch.Tensor, torch.Tensor],
        dataloader_idx: int,
    ) -> dict[str, torch.Tensor | None]:
        """Format the datapipe output with registry keys for scvi-tools."""
        X, obs = batch

        X_KEY: str = "X"
        BATCH_KEY: str = "batch"
        LABELS_KEY: str = "labels"

        return {
            X_KEY: X,
            BATCH_KEY: obs,
            LABELS_KEY: None,
        }
