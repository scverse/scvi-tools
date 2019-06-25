import pickle
import os
import logging
import numpy as np

from scvi.dataset.dataset import GeneExpressionDataset


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


class SyntheticRandomDataset(GeneExpressionDataset):
    # The exact random parameters taken by romain with simlr labels
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
        logging.info("mu " + str(mu) + " theta " + str(theta) + " r " + str(r) + " p " + str(
            p) + " dropout " + str(dropout))

        self.save_path = save_path
        self.url = 'https://github.com/YosefLab/scVI-data/raw/master/random_metadata.pickle'
        self.download_name = 'random_metadata.pickle'
        self.download()

        self.simlr_metadata = pickle.load(
            open(os.path.join(self.save_path, 'random_metadata.pickle'), 'rb'))
        labels_simlr = self.simlr_metadata['clusters']

        super().__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(X_train, labels=labels_simlr)
        )


class SyntheticDatasetCorr(GeneExpressionDataset):
    def __init__(self, n_cells_cluster=200, n_clusters=3,
                 n_genes_high=25, n_overlap=0,
                 weight_high=4e-2, weight_low=1e-2,
                 lam_0=50., n_batches=1):
        """
        We define technical zeros of the synthetic dataset as the zeros that result from
        highly expressed genes (relatively to the considered cell) and the biological zeros as the
        rest of the zeros
        :param n_cells_cluster: Number of cells in each cluster
        :param n_clusters: Number of cell clusters
        :param n_genes_high: Number of highly expressed genes for each cluster.
        Other genes are still expressed but at a lower expression
        The level expression for highly expressed genes follow the same distribution
        :param n_genes_total: Number of genes in total
        :param weight_high: Level weight for highly expressed genes
        :param weight_low: Level weight for lowly expressed genes
        :param lam_0: Proportionality in the Poisson distributions parameter
        :param n_batches: (useless for now) Number of batches in dataset
        """
        assert n_batches == 1
        np.random.seed(0)

        if n_overlap == 0:
            n_genes_total = n_clusters * n_genes_high
        else:
            n_genes_total = n_clusters * (n_genes_high - n_overlap) + n_overlap

        if n_genes_total % n_clusters > 0:
            logging.info("Warning, clusters have inequal sizes")

        if n_genes_high > (n_genes_total // n_clusters):
            logging.info("Overlap of", n_genes_high - (n_genes_total // n_clusters), "genes")

        # Generate data before dropout
        batch_size = n_cells_cluster * n_clusters
        self.n_cells_cluster = n_cells_cluster
        self.exprs_param = np.ones((n_batches, batch_size, n_genes_total))
        self.n_clusters = n_clusters
        self.batch_size = batch_size
        self.n_genes_total = n_genes_total
        self.n_genes_high = n_genes_high
        self.is_highly_exp = np.zeros((self.n_batches, self.batch_size, self.n_genes_total),
                                      dtype=np.bool)
        self.probas_zero_bio_tech_high = None
        self.probas_zero_bio_tech_low = None
        labels = np.ones((n_batches, batch_size, 1))

        # For each cell cluster, some genes have a high expression, the rest
        # has a low expression. The scope of high expression genes "moves"
        # with the cluster
        for cluster in range(n_clusters):
            labels[:, cluster * n_cells_cluster:(cluster + 1) * n_cells_cluster, :] = cluster

            ind_first_gene_cluster = cluster * (n_genes_high - n_overlap)
            ind_last_high_gene_cluster = ind_first_gene_cluster + n_genes_high
            self.is_highly_exp[:,
                               cluster * n_cells_cluster:(cluster + 1) * n_cells_cluster,
                               ind_first_gene_cluster:ind_last_high_gene_cluster] = True
            # Weights in a cluster to create highly-expressed and low-expressed genes
            weights = weight_low * np.ones((n_genes_total,))
            weights[ind_first_gene_cluster:ind_last_high_gene_cluster] = weight_high
            weights /= weights.sum()

            self.exprs_param[:,
                             cluster*n_cells_cluster:(cluster+1)*n_cells_cluster,
                             :] = lam_0 * weights

        logging.info('Poisson Params extremal values: ', self.exprs_param.min(), self.exprs_param.max())
        expression_mat = np.random.poisson(self.exprs_param)
        self.poisson_zero = np.exp(-self.exprs_param)
        new_data = self.mask(expression_mat)
        logging.info((new_data == 0).sum())
        super().__init__(
            *GeneExpressionDataset.get_attributes_from_list(new_data, list_labels=labels),
            gene_names=np.arange(n_genes_total).astype(np.str))

    def mask(self, data):
        return data


class ZISyntheticDatasetCorr(SyntheticDatasetCorr):
    def __init__(self, dropout_coef_high=0.05, lam_dropout_high=0.,
                 dropout_coef_low=0.08, lam_dropout_low=0., **kwargs):
        assert max(dropout_coef_high, dropout_coef_low) < 1
        self.dropout_coef_high = dropout_coef_high
        self.lam_dropout_high = lam_dropout_high
        self.dropout_coef_low = dropout_coef_low
        self.lam_dropout_low = lam_dropout_low
        self.p_dropout = None
        self.zi_zero = None
        self.is_technical = None
        super(ZISyntheticDatasetCorr, self).__init__(**kwargs)

    def mask(self, data):
        self.p_dropout = self.dropout_coef_low * np.exp(
            -self.lam_dropout_low * (self.exprs_param ** 2))

        self.p_dropout[self.is_highly_exp] = (
            self.dropout_coef_high
            * np.exp(
                -self.lam_dropout_high * (self.exprs_param ** 2)
              ))[self.is_highly_exp]

        # Probability of failure
        mask = np.random.binomial(n=1, p=1 - self.p_dropout,
                                  size=(self.n_batches, self.batch_size, self.n_genes_total))
        self.is_technical = mask == 0
        self.zi_zero = self.p_dropout

        # Probas of ZI = or != 0 and NB = or != 0 for highly and lowly expressed genes
        self.probas_zero_bio_tech_high = [
            [np.mean(((1 - self.poisson_zero) * (1 - self.zi_zero))[self.is_highly_exp]),
             # ZI != 0, NB != 0
             np.mean(((1 - self.poisson_zero) * self.zi_zero)[self.is_highly_exp])],
            # ZI = 0, NB != 0
            [np.mean((self.poisson_zero * (1 - self.zi_zero))[self.is_highly_exp]),
             # ZI != 0, NB = 0
             np.mean((self.poisson_zero * self.zi_zero)[self.is_highly_exp])]]  # ZI = 0, NB = 0
        self.probas_zero_bio_tech_high = np.array(self.probas_zero_bio_tech_high).reshape((2, 2))

        self.probas_zero_bio_tech_low = [
            [np.mean(((1 - self.poisson_zero) * (1 - self.zi_zero))[~self.is_highly_exp]),
             # ZI != 0, NB != 0
             np.mean(((1 - self.poisson_zero) * self.zi_zero)[~self.is_highly_exp])],
            # ZI = 0, NB != 0
            [np.mean((self.poisson_zero * (1 - self.zi_zero))[~self.is_highly_exp]),
             # ZI != 0, NB = 0
             np.mean((self.poisson_zero * self.zi_zero)[~self.is_highly_exp])]]  # ZI = 0, NB = 0
        self.probas_zero_bio_tech_low = np.array(self.probas_zero_bio_tech_low).reshape((2, 2))
        assert np.abs(self.probas_zero_bio_tech_high.sum() - 1) <= 1e-8
        assert np.abs(self.probas_zero_bio_tech_low.sum() - 1) <= 1e-8
        return data * mask.astype(np.float32)
