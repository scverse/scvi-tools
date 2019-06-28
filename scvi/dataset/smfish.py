import logging
import os

import loompy
import numpy as np

from scvi.dataset.dataset import DownloadableDataset

logger = logging.getLogger(__name__)


class SmfishDataset(DownloadableDataset):
    def __init__(self, save_path="data/", cell_type_level="major", preprocess=True):
        self._preprocess_data = preprocess
        self.cell_type_level = cell_type_level
        super().__init__(
            "http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom",
            "osmFISH_SScortex_mouse_all_cell.loom",
            save_path,
            delayed_populating=False,
        )

    def populate(self):
        data = self.load_from_disk()
        self.populate_from_data(**data)
        if self._preprocess_data:
            self.preprocess()
        self.filter_cells_by_count()

    def load_from_disk(self):
        logger.info("Loading smFISH dataset")
        ds = loompy.connect(os.path.join(self.save_path, self.filenames[0]))
        gene_names = np.char.upper(ds.ra["Gene"].astype(np.str))

        labels = ds.ca["ClusterID"].reshape(-1, 1)
        tmp_cell_types = np.array(ds.ca["ClusterName"])

        u_labels, u_index = np.unique(labels.ravel(), return_index=True)
        cell_types = ["" for _ in range(max(u_labels) + 1)]
        for i, index in zip(u_labels, u_index):
            cell_types[i] = tmp_cell_types[index]
        cell_types = np.array(cell_types, dtype=np.str)

        x_coord, y_coord = ds.ca["X"], ds.ca["Y"]
        x_coord = x_coord.reshape((-1, 1))
        y_coord = y_coord.reshape((-1, 1))
        data = ds[:, :].T

        return {
            "X": data,
            "labels": labels,
            "gene_names": gene_names,
            "cell_types": cell_types,
            "cell_attributes_dict": {"x_coord": x_coord, "y_coord": y_coord},
        }

    def preprocess(self):
        major_clusters = dict(
            [
                ((3, 2), "Astrocytes"),
                ((7, 26), "Endothelials"),
                ((18, 17, 14, 19, 15, 16, 20), "Inhibitory"),
                ((29, 28), "Microglias"),
                ((32, 33, 30, 22, 21), "Oligodendrocytes"),
                ((9, 8, 10, 6, 5, 4, 12, 1, 13), "Pyramidals"),
            ]
        )
        if self.cell_type_level == "major":
            self.map_cell_types(major_clusters)
            self.filter_cell_types(
                [
                    "Astrocytes",
                    "Endothelials",
                    "Inhibitory",
                    "Microglias",
                    "Oligodendrocytes",
                    "Pyramidals",
                ]
            )
        self.remap_categorical_attributes()
