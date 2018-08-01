import numpy as np
import loompy
from scvi.dataset import GeneExpressionDataset


class SmfishDatasetCustom(GeneExpressionDataset):
    def __init__(self, save_path='data/'):

        self.download_name = 'osmFISH_SScortex_mouse_all_cell.loom'
        self.save_path = save_path
        self.url = 'http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom'

        data, labels, cell_types, x_coord, y_coord, gene_names = self.download_and_preprocess()
        cell_types = ["astrocytes_ependymal",
                      "endothelial-mural",
                      "interneurons",
                      "microglia",
                      "oligodendrocytes",
                      "pyramidal"]
        super(SmfishDatasetCustom, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                data,
                labels=labels), gene_names=gene_names,
            x_coord=x_coord, y_coord=y_coord, cell_types=cell_types)

    def preprocess(self):
        print("Preprocessing smFISH dataset")
        ds = loompy.connect(self.save_path + self.download_name)
        # select = ds[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene
        gene_names = ds.ra['Gene']
        # print(gene_names)
        labels = ds.ca['ClusterID']
        to_keep = []
        to_discard = 0
        for n_label in range(len(labels)):
            if labels[n_label] not in [0, 27, 26, 24, 31]:
                to_keep.append(n_label)
            else:
                to_discard += 1
        select = np.array(to_keep)
        labels, cell_types = np.array(labels), np.array(ds.ca['ClusterName'])
        labels = np.reshape(labels, (labels.shape[0], 1))[select]
        # cell_types = np.reshape(cell_types, (cell_types.shape[0], 1))[select]
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
        new_labels = []
        new_cell_types = []
        for label in labels:
            if label in major_clusters_fish['Astrocytes']:
                new_labels.append(0)
                new_cell_types.append('astrocytes')
            if label in major_clusters_fish['Endothelial']:
                new_labels.append(1)
                new_cell_types.append('endothelial')
            if label in major_clusters_fish['Inhibitory']:
                new_labels.append(2)
                new_cell_types.append('inhibitory')
            if label in major_clusters_fish['Microglia']:
                new_labels.append(3)
                new_cell_types.append('microglia')
            if label in major_clusters_fish['Oligodendrocytes']:
                new_labels.append(4)
                new_cell_types.append('oligodendrocytes')
            if label in major_clusters_fish['Excitatory']:
                new_labels.append(5)
                new_cell_types.append('excitatory')

        labels, cell_types = np.array(new_labels).reshape(-1, 1), np.array(new_cell_types).reshape(-1, 1)

        x_coord, y_coord = np.array(ds.ca['X']), np.array(ds.ca['Y'])
        x_coord = np.reshape(x_coord, (x_coord.shape[0], 1))[select]
        y_coord = np.reshape(y_coord, (y_coord.shape[0], 1))[select]

        data = ds[:, select].T  # change matrix to cells by genes

        select = data.sum(axis=1) > 0  # Take out cells that doesn't express any gene
        data = data[select, :]
        labels = labels[select]
        cell_types = cell_types[select]
        x_coord = x_coord[select]
        y_coord = y_coord[select]

        ds.close()

        print("Finished preprocessing smFISH dataset")
        return data, labels, cell_types, x_coord, y_coord, gene_names
