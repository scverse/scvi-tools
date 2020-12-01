from typing import Optional

import numpy as np
from torch.utils.data import DataLoader
from ._scvi_dataloader import ScviDataLoader
from itertools import cycle


class ConcatDataLoader(DataLoader):
    def __init__(
        self,
        adata,
        indices_list,
        shuffle=False,
        batch_size=128,
        data_and_attributes: Optional[dict] = None,
        **data_loader_kwargs,
    ):
        self.dataloaders = []
        for indices in indices_list:
            self.dataloaders.append(
                ScviDataLoader(
                    adata,
                    indices=indices,
                    shuffle=True,
                    batch_size=batch_size,
                    data_and_attributes=data_and_attributes,
                    **data_loader_kwargs,
                )
            )
        lens = [len(dl) for dl in self.dataloaders]
        self.largest_dl = self.dataloaders[np.argmax(lens)]
        super().__init__(self.largest_dl, **data_loader_kwargs)

    def __len__(self):
        return len(self.largest_dl)

    def __iter__(self):
        iter_list = [
            cycle(dl) if dl != self.largest_dl else dl for dl in self.dataloaders
        ]
        return zip(*iter_list)
