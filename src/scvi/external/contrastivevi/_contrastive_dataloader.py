import warnings
from itertools import cycle

from scvi import settings
from scvi.data import AnnDataManager
from scvi.dataloaders import ConcatDataLoader


class _ContrastiveIterator:
    """Iterator for background and target dataloader pair in contrastive analysis.

    Each iteration of this iterator returns a dictionary with two elements:
    "background", containing one batch of data from the background dataloader, and
    "target", containing one batch of data from the target dataloader.
    """

    def __init__(self, background, target):
        self.background = iter(background)
        self.target = iter(target)

    def __iter__(self):
        return self

    def __next__(self):
        bg_samples = next(self.background)
        tg_samples = next(self.target)
        return {"background": bg_samples, "target": tg_samples}


class ContrastiveDataLoader(ConcatDataLoader):
    """Dataloader to load background and target data for contrastive analysis.

    Each iteration of the dataloader returns a dictionary containing background and
    target data points, indexed by "background" and "target", respectively.

    Parameters
    ----------
        adata_manager
            :class:`~scvi.data.AnnDataManager` object that has been created via
            ``setup_anndata``.
        background_indices
            Indices for background samples in the adata.
        target_indices
            Indices for target samples in the adata.
        shuffle
            Whether the data should be shuffled.
        batch_size
            Mini-batch size to load for background and target data.
        data_and_attributes
            Dictionary with keys representing keys in data registry
            (``adata_manager.data_registry``) and value equal to desired numpy loading
            type (later made into torch tensor). If ``None``, defaults to all registered
            data.
        drop_last
            If ``int``, drops the last batch if its length is less than
            ``drop_last``. If ``drop_last == True``, drops last non-full batch.
            If ``drop_last == False``, iterate over all batches.
        distributed_sampler
            ``EXPERIMENTAL`` Whether to use :class:`~scvi.dataloaders.BatchDistributedSampler` as
            the sampler. If `True`, `sampler` must be `None`. Not applicable here.
        load_sparse_tensor
            ``EXPERIMENTAL`` If ``True``, loads data with sparse CSR or CSC layout as a
            :class:`~torch.Tensor` with the same layout. Can lead to speedups in data
            transfers to GPUs, depending on the sparsity of the data. Not applicable
            here.
        **data_loader_kwargs: Keyword arguments for `torch.utils.data.DataLoader`.
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        background_indices: list[int],
        target_indices: list[int],
        shuffle: bool = False,
        batch_size: int = 128,
        data_and_attributes: dict | None = None,
        drop_last: bool | int = False,
        distributed_sampler: bool = False,
        load_sparse_tensor: bool = False,
        **data_loader_kwargs,
    ) -> None:
        super().__init__(
            adata_manager=adata_manager,
            indices_list=[background_indices, target_indices],
            shuffle=shuffle,
            batch_size=batch_size,
            data_and_attributes=data_and_attributes,
            drop_last=drop_last,
            **data_loader_kwargs,
        )
        self.background_indices = background_indices
        self.target_indices = target_indices
        if distributed_sampler:
            warnings.warn(
                (
                    "distributed_sampler=True is not implemented for "
                    "ContrastiveDataLoader. Setting distributed_sampler=False"
                ),
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )
            distributed_sampler = False
        if load_sparse_tensor:
            warnings.warn(
                (
                    "load_sparse_tensor=True is not implemented for "
                    "ContrastiveDataLoader. Setting load_sparse_tensor=False"
                ),
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )
            load_sparse_tensor = False
        self.distributed_sampler = distributed_sampler
        self.load_sparse_tensor = load_sparse_tensor

    def __iter__(self):
        """Iter method for contrastive dataloader.

        Will iter over the dataloader with the most data while cycling through
        the data in the other dataloader.
        """
        iter_list = [cycle(dl) if dl != self.largest_dl else dl for dl in self.dataloaders]
        return _ContrastiveIterator(background=iter_list[0], target=iter_list[1])
