import logging
import os

import pandas as pd

from scvi.dataset.dataset import DownloadableDataset

logger = logging.getLogger(__name__)


class SeqfishDataset(DownloadableDataset):
    def __init__(self, save_path: str = "data/", delayed_populating: bool = False):
        super().__init__(
            urls="https://www.cell.com/cms/attachment/2080562255/2072099886/mmc6.xlsx",
            filenames="SeqFISH.xlsx",
            save_path=save_path,
            delayed_populating=delayed_populating,
        )

    def populate(self):
        logger.info("Preprocessing dataset")

        xl = pd.ExcelFile(os.path.join(self.save_path, self.filenames[0]))
        ds = xl.parse("Hippocampus Counts")  # They also used cell by genes

        logger.info("Finished preprocessing dataset")
        self.populate_from_data(ds.values[:, 1:].astype(int))
