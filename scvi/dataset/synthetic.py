import pickle
import os
import numpy as np

from . import GeneExpressionDataset


class SyntheticDataset(GeneExpressionDataset):
    def __init__(self, batch_size=200, nb_genes=100, n_batches=2, n_labels=3):
        # Generating samples according to a ZINB process
        data = np.random.negative_binomial(5, 0.3, size=(n_batches, batch_size, nb_genes))
        mask = np.random.binomial(n=1, p=0.7, size=(n_batches, batch_size, nb_genes))
        newdata = (data * mask)  # We put the batch index first
        labels = np.random.randint(0, n_labels, size=(n_batches, batch_size, 1))
        super().__init__(
            *GeneExpressionDataset.get_attributes_from_list(newdata, list_labels=labels),
            gene_names=np.arange(nb_genes).astype(np.str)
        )


class SyntheticRandomDataset(GeneExpressionDataset):  # The exact random parameters taken by romain with simlr labels
    def __init__(self, mu=4.0, theta=2.0, dropout=0.7, save_path='data/'):
        np.random.seed(0)  # simlr attributes computed with this seed.
        p = mu / (mu + theta)
        r = theta
        shape_train = (2000, 10)
        l_train = np.random.gamma(r, p / (1 - p), size=shape_train)
        X_train = np.random.poisson(l_train)
        X_train *= np.random.binomial(1, 1 - dropout, size=shape_train)
        indices_to_keep = (X_train.sum(axis=1) > 0).ravel()
        X_train = X_train[indices_to_keep]
        print("mu " + str(mu) + " theta " + str(theta) + " r " + str(r) + " p " + str(p) + " dropout " + str(dropout))

        self.save_path = save_path
        self.url = 'https://github.com/YosefLab/scVI-data/raw/master/random_metadata.pickle'
        self.download_name = 'random_metadata.pickle'
        self.download()

        self.simlr_metadata = pickle.load(open(os.path.join(self.save_path, 'random_metadata.pickle'), 'rb'))
        labels_simlr = self.simlr_metadata['clusters']

        super().__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(X_train, labels=labels_simlr)
        )


class SyntheticDatasetCorr(GeneExpressionDataset):
    def __init__(self, n_cells_cluster=100, n_clusters=5,
                 n_genes_high=15, n_genes_total=50,
                 weight_high=0.1, weight_low=0.01,
                 lam_0=100., n_batches=1, n_labels=1, mode='mixed'):
        """

        :param n_cells_cluster: Number of cells in each cluster
        :param n_clusters: Number of cell clusters
        :param n_genes_high: Number of highly expressed genes for each cluster.
        Other genes are still expressed but at a lower expression
        The level expression for highly expressed genes follow the same distribution

        :param n_genes_total: Number of genes in total
        :param weight_high: Level weight for highly expressed genes
        :param weight_low: Level weight for lowly expressed genes
        :param lam_0: Proportionality in the Poisson distributions parameter
        :param p_dropout: Probability of dropout (only makes sense if mode='zi'
        or mode='mixed'
        :param n_batches: (useless for now) Number of batches in dataset
        :param n_labels: (useless for now) Number of labels
        :param mode: Should be either
            - 'nb': mode without zero inflation
            - 'zi': mode including zero inflation for ALL genes
            -'mixed': mode in which SOME genes are Zero-Inflated while some other are not
            The ratio_genes_zi parameter allow to control the proportion of genes that are going
            to be ZI
        :param ratio_genes_zi: Proportion of ZI genes
        """

        assert mode in ['mixed', 'zi', 'nb'], 'Mode {} not recognized'.format(mode)
        np.random.seed(0)

        if n_genes_total % n_clusters > 0:
            print("Warning, clusters have inequal sizes")

        if n_genes_high > (n_genes_total // n_clusters):
            print("Overlap of", n_genes_high - (n_genes_total // n_clusters), "genes")

        # Generate data before dropout
        batch_size = n_cells_cluster * n_clusters
        self.exprs_param = np.ones((n_batches, batch_size, n_genes_total))

        self.batch_size = batch_size
        self.n_batches = n_batches
        self.n_genes_total = n_genes_total

        labels = np.ones((n_batches, batch_size, 1))

        # For each cell cluster, some genes have a high expression, the rest
        # has a low expression. The scope of high expression genes "moves"
        # with the cluster
        for cluster in range(n_clusters):

            labels[:, cluster * n_cells_cluster:(cluster + 1) * n_cells_cluster, :] = cluster


            ind_first_gene_cluster = cluster * (n_genes_total // n_clusters)
            ind_last_high_gene_cluster = ind_first_gene_cluster + n_genes_high

            # Weights in a cluster to create highly-expressed and low-expressed genes
            weights = weight_low * np.ones((n_genes_total,))
            weights[ind_first_gene_cluster:ind_last_high_gene_cluster] = weight_high
            weights /= weights.sum()

            self.exprs_param[:, cluster * n_cells_cluster:(cluster + 1) * n_cells_cluster, :] = lam_0 * weights

        # Apply dropout depending on the mode
        expression_mat = np.random.poisson(self.exprs_param)

        new_data = self.mask(expression_mat)
        print((new_data== 0).sum())
        super().__init__(
            *GeneExpressionDataset.get_attributes_from_list(new_data, list_labels=labels),
            gene_names=np.arange(n_genes_total).astype(np.str))

    def mask(self, data):
        return data


class ZISyntheticDatasetCorr(SyntheticDatasetCorr):
    def __init__(self, dropout_coef=0.5, lam_dropout=1.0, **kwargs):
        assert dropout_coef < 1
        self.dropout_coef = dropout_coef
        self.lam_dropout = lam_dropout
        self.is_technical = None
        super(ZISyntheticDatasetCorr, self).__init__(**kwargs)

    def mask(self, data):
        p_dropout = self.dropout_coef * np.exp(-self.lam_dropout * (self.exprs_param**2))
        # Probability of failure
        mask = np.random.binomial(n=1, p=1 - p_dropout,
                                  size=(self.n_batches, self.batch_size, self.n_genes_total))
        self.is_technical = (mask == 0)
        return data * mask.astype(np.float32)
