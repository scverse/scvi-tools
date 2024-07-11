from collections.abc import Sequence

from anndata import AnnData


def get_aggregated_posterior(
    self,
    adata: AnnData | None = None,
    indices: Sequence[int] | None = None,
    batch_size: int | None = None,
):
    adata = self._validate_anndata(adata)  # need to maybe inherit form vaemixin for now?
    dataloader = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

    for _tensors in dataloader:
        pass


# def differential_abundance():
