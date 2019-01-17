import numpy as np
import loompy
from .dataset import GeneExpressionDataset
import os


class SmfishDataset(GeneExpressionDataset):
    def __init__(self, save_path='data/', cell_type_level="major"):

        self.download_name = 'osmFISH_SScortex_mouse_all_cell.loom'
        self.save_path = save_path
        self.url = 'http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom'
        self.cell_type_level = cell_type_level
        data, labels, gene_names, cell_types, x_coord, y_coord = self.download_and_preprocess()
        super().__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                data,
                labels=labels), gene_names=gene_names,
            x_coord=x_coord, y_coord=y_coord)

    def preprocess(self):
        print("Preprocessing smFISH dataset")
        ds = loompy.connect(os.path.join(self.save_path, self.download_name))
        gene_names = ds.ra['Gene']
        if self.cell_type_level == "minor":
            select = ds[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene

            labels, cell_types = np.array(ds.ca['ClusterID']), np.array(ds.ca['ClusterName'])
            labels = np.reshape(labels, (labels.shape[0], 1))[select]
            cell_types = np.reshape(cell_types, (cell_types.shape[0], 1))[select]

        elif self.cell_type_level == "major":
            major_clusters_fish = {
                'Inhibitory': [18, 17, 14, 19, 15, 16, 20],
                'Excitatory': [9, 8, 10, 6, 5, 4, 12, 1, 13],
                'Astrocytes': [3, 2],
                'Oligodendrocytes': [32, 33, 30, 22, 21],
                'Microglia': [29, 28],
                'Choroid plexus': [24],
                'Ependimal': [27],
                'Pericytes': [31],
                'Endothelial': [7, 25],
                'VSM': [25]
            }
            labels = ds.ca['ClusterID']
            to_keep = []
            new_labels = []
            for n_label in range(len(labels)):
                if labels[n_label] not in [0, 27, 26, 24, 31]:
                    to_keep.append(n_label)
                if labels[n_label] in major_clusters_fish['Astrocytes']:
                    new_labels.append(0)
                elif labels[n_label] in major_clusters_fish['Endothelial']:
                    new_labels.append(1)
                elif labels[n_label] in major_clusters_fish['Inhibitory']:
                    new_labels.append(2)
                elif labels[n_label] in major_clusters_fish['Microglia']:
                    new_labels.append(3)
                elif labels[n_label] in major_clusters_fish['Oligodendrocytes']:
                    new_labels.append(4)
                elif labels[n_label] in major_clusters_fish['Excitatory']:
                    new_labels.append(5)

            select = np.array(to_keep)
            labels, cell_types = np.array(new_labels), np.array(ds.ca['ClusterName'])
            labels = np.reshape(labels, (labels.shape[0], 1))
            cell_types = np.reshape(cell_types, (cell_types.shape[0], 1))[select]

        x_coord, y_coord = np.array(ds.ca['X']), np.array(ds.ca['Y'])
        x_coord = np.reshape(x_coord, (x_coord.shape[0], 1))[select]
        y_coord = np.reshape(y_coord, (y_coord.shape[0], 1))[select]

        data = ds[:, select].T

        print("Finished preprocessing smFISH dataset")
        return data, labels, gene_names, cell_types, x_coord, y_coord
