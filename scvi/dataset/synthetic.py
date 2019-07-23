import logging
import os
import pickle

import numpy as np

from scvi.dataset.dataset import GeneExpressionDataset, DownloadableDataset

logger = logging.getLogger(__name__)


class SyntheticDataset(GeneExpressionDataset):
    def __init__(
        self,
        batch_size: int = 200,
        nb_genes: int = 100,
        n_batches: int = 2,
        n_labels: int = 3,
    ):
        super().__init__()
        # Generating samples according to a ZINB process
        data = np.random.negative_binomial(
            5, 0.3, size=(n_batches, batch_size, nb_genes)
        )
        mask = np.random.binomial(n=1, p=0.7, size=(n_batches, batch_size, nb_genes))
        data = data * mask  # We put the batch index first
        labels = np.random.randint(0, n_labels, size=(n_batches, batch_size, 1))
        cell_types = ["undefined_%d" % i for i in range(n_labels)]

        self.populate_from_per_batch_list(
            data,
            labels_per_batch=labels,
            gene_names=np.arange(nb_genes).astype(np.str),
            cell_types=cell_types,
        )
        # clear potentially unused cell_types
        self.remap_categorical_attributes()


class SyntheticRandomDataset(DownloadableDataset):
    FILENAME = "random_metadata.pickle"

    # The exact random parameters taken by romain with simlr labels
    def __init__(
        self,
        mu: float = 4.0,
        theta: float = 2.0,
        dropout: float = 0.7,
        save_path: str = "data/",
    ):
        self.mu = mu
        self.theta = theta
        self.dropout = dropout
        self.simlr_metadata = None
        super().__init__(
            urls="https://github.com/YosefLab/scVI-data/raw/master/random_metadata.pickle",
            filenames=SyntheticRandomDataset.FILENAME,
            save_path=save_path,
        )

    def populate(self):
        np.random.seed(0)  # simlr attributes computed with this seed.
        p = self.mu / (self.mu + self.theta)
        r = self.theta
        shape_train = (2000, 10)
        l_train = np.random.gamma(r, p / (1 - p), size=shape_train)
        X_train = np.random.poisson(l_train)
        X_train *= np.random.binomial(1, 1 - self.dropout, size=shape_train)
        indices_to_keep = (X_train.sum(axis=1) > 0).ravel()
        X_train = X_train[indices_to_keep]
        logger.info(
            "mu "
            + str(self.mu)
            + " theta "
            + str(self.theta)
            + " r "
            + str(r)
            + " p "
            + str(p)
            + " dropout "
            + str(self.dropout)
        )

        self.simlr_metadata = pickle.load(
            open(os.path.join(self.save_path, SyntheticRandomDataset.FILENAME), "rb")
        )
        labels_simlr = self.simlr_metadata["clusters"]

        self.populate_from_data(X_train, labels=labels_simlr)


class SyntheticDatasetCorr(GeneExpressionDataset):
    def __init__(
        self,
        n_cells_cluster: int = 200,
        n_clusters: int = 3,
        n_genes_high: int = 25,
        n_overlap: int = 0,
        weight_high: float = 4e-2,
        weight_low: float = 1e-2,
        lam_0: float = 50.0,
    ):
        """
        We define technical zeros of the synthetic dataset as the zeros that result from
        highly expressed genes (relatively to the considered cell) and the biological zeros as the
        rest of the zeros
        :param n_cells_cluster: Number of cells in each cluster
        :param n_clusters: Number of cell clusters
        :param n_genes_high: Number of highly expressed genes for each cluster.
        Other genes are still expressed but at a lower expression
        The level expression for highly expressed genes follow the same distribution
        :param weight_high: Level weight for highly expressed genes
        :param weight_low: Level weight for lowly expressed genes
        :param lam_0: Proportionality in the Poisson distributions parameter
        """
        super().__init__()
        n_batches = 1
        np.random.seed(0)

        if n_overlap == 0:
            n_genes_total = n_clusters * n_genes_high
        else:
            n_genes_total = n_clusters * (n_genes_high - n_overlap) + n_overlap

        if n_genes_total % n_clusters > 0:
            logger.info("Warning, clusters have inequal sizes")

        if n_genes_high > (n_genes_total // n_clusters):
            logger.info(
                "Overlap of {} genes".format(
                    n_genes_high - (n_genes_total // n_clusters)
                )
            )

        # Generate data before dropout
        batch_size = n_cells_cluster * n_clusters
        self.n_cells_cluster = n_cells_cluster
        self.exprs_param = np.ones((n_batches, batch_size, n_genes_total))
        self.n_clusters = n_clusters
        self.batch_size = batch_size
        self.n_genes_total = n_genes_total
        self.n_genes_high = n_genes_high
        self.is_highly_exp = np.zeros(
            (n_batches, self.batch_size, self.n_genes_total), dtype=np.bool
        )
        self.probas_zero_bio_tech_high = None
        self.probas_zero_bio_tech_low = None
        labels = np.ones((n_batches, batch_size, 1))

        # For each cell cluster, some genes have a high expression, the rest
        # has a low expression. The scope of high expression genes "moves"
        # with the cluster
        for cluster in range(n_clusters):
            labels[
                :, cluster * n_cells_cluster : (cluster + 1) * n_cells_cluster, :
            ] = cluster

            ind_first_gene_cluster = cluster * (n_genes_high - n_overlap)
            ind_last_high_gene_cluster = ind_first_gene_cluster + n_genes_high
            self.is_highly_exp[
                :,
                cluster * n_cells_cluster : (cluster + 1) * n_cells_cluster,
                ind_first_gene_cluster:ind_last_high_gene_cluster,
            ] = True
            # Weights in a cluster to create highly-expressed and low-expressed genes
            weights = weight_low * np.ones((n_genes_total,))
            weights[ind_first_gene_cluster:ind_last_high_gene_cluster] = weight_high
            weights /= weights.sum()

            self.exprs_param[
                :, cluster * n_cells_cluster : (cluster + 1) * n_cells_cluster, :
            ] = (lam_0 * weights)

        logger.info(
            "Poisson Params extremal values: {min}, {max}".format(
                min=self.exprs_param.min(), max=self.exprs_param.max()
            )
        )
        expression_mat = np.random.poisson(self.exprs_param)
        self.poisson_zero = np.exp(-self.exprs_param)
        new_data = self.mask(expression_mat)
        logger.info((new_data == 0).sum())

        self.populate_from_per_batch_list(
            new_data,
            labels_per_batch=labels,
            gene_names=np.arange(n_genes_total).astype(np.str),
        )

    def mask(self, data):
        return data


class ZISyntheticDatasetCorr(SyntheticDatasetCorr):
    def __init__(
        self,
        dropout_coef_high: float = 0.05,
        lam_dropout_high: float = 0.0,
        dropout_coef_low: float = 0.08,
        lam_dropout_low: float = 0.0,
        n_cells_cluster: int = 200,
        n_clusters: int = 3,
        n_genes_high: int = 25,
        n_overlap: int = 0,
        weight_high: float = 4e-2,
        weight_low: float = 1e-2,
        lam_0: float = 50.0,
    ):
        assert max(dropout_coef_high, dropout_coef_low) < 1
        self.dropout_coef_high = dropout_coef_high
        self.lam_dropout_high = lam_dropout_high
        self.dropout_coef_low = dropout_coef_low
        self.lam_dropout_low = lam_dropout_low
        self.p_dropout = None
        self.zi_zero = None
        self.is_technical = None
        super().__init__(
            n_cells_cluster=n_cells_cluster,
            n_clusters=n_clusters,
            n_genes_high=n_genes_high,
            n_overlap=n_overlap,
            weight_high=weight_high,
            weight_low=weight_low,
            lam_0=lam_0,
        )

    def mask(self, data):
        self.p_dropout = self.dropout_coef_low * np.exp(
            -self.lam_dropout_low * (self.exprs_param ** 2)
        )

        self.p_dropout[self.is_highly_exp] = (
            self.dropout_coef_high
            * np.exp(-self.lam_dropout_high * (self.exprs_param ** 2))
        )[self.is_highly_exp]

        # Probability of failure
        mask = np.random.binomial(n=1, p=1 - self.p_dropout, size=data.shape)
        self.is_technical = mask == 0
        self.zi_zero = self.p_dropout

        # Probas of ZI = or != 0 and NB = or != 0 for highly and lowly expressed genes
        self.probas_zero_bio_tech_high = [
            [
                np.mean(
                    ((1 - self.poisson_zero) * (1 - self.zi_zero))[self.is_highly_exp]
                ),
                # ZI != 0, NB != 0
                np.mean(((1 - self.poisson_zero) * self.zi_zero)[self.is_highly_exp]),
            ],
            # ZI = 0, NB != 0
            [
                np.mean((self.poisson_zero * (1 - self.zi_zero))[self.is_highly_exp]),
                # ZI != 0, NB = 0
                np.mean((self.poisson_zero * self.zi_zero)[self.is_highly_exp]),
            ],
        ]  # ZI = 0, NB = 0
        self.probas_zero_bio_tech_high = np.asarray(
            self.probas_zero_bio_tech_high
        ).reshape((2, 2))

        self.probas_zero_bio_tech_low = [
            [
                np.mean(
                    ((1 - self.poisson_zero) * (1 - self.zi_zero))[~self.is_highly_exp]
                ),
                # ZI != 0, NB != 0
                np.mean(((1 - self.poisson_zero) * self.zi_zero)[~self.is_highly_exp]),
            ],
            # ZI = 0, NB != 0
            [
                np.mean((self.poisson_zero * (1 - self.zi_zero))[~self.is_highly_exp]),
                # ZI != 0, NB = 0
                np.mean((self.poisson_zero * self.zi_zero)[~self.is_highly_exp]),
            ],
        ]  # ZI = 0, NB = 0
        self.probas_zero_bio_tech_low = np.asarray(
            self.probas_zero_bio_tech_low
        ).reshape((2, 2))
        assert np.abs(self.probas_zero_bio_tech_high.sum() - 1) <= 1e-8
        assert np.abs(self.probas_zero_bio_tech_low.sum() - 1) <= 1e-8
        return data * mask.astype(np.float32)
