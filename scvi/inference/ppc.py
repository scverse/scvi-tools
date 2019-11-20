from typing import Dict
import numpy as np
from scvi.inference import Posterior
from sklearn.decomposition import FactorAnalysis, PCA
from sklearn.preprocessing import StandardScaler
from scipy import linalg
import pandas as pd
from tqdm.auto import tqdm


class PosteriorPredictiveCheck:
    """Posterior predictive checks for comparing scVI models
    and ability to compare to any other models
    """

    def __init__(
        self, posteriors_dict: Dict[str, Posterior], n_samples: int = 10, batch_size=32
    ):
        """
        Args:
            posteriors_dict (Dict[str, Posterior]): dictionary of Posterior objects fit on the same dataset
            n_samples (int, optional): Number of posterior predictive samples. Defaults to 10.
        """
        self.posteriors = posteriors_dict
        self.dataset = posteriors_dict[next(iter(posteriors_dict.keys()))].gene_dataset
        self.raw_counts = None
        self.posterior_predictive_samples = {}
        self.n_samples = n_samples
        self.models = {}
        self.metrics = {}
        self.raw_metrics = {}
        self.batch_size = batch_size

        self.store_posterior_samples()

    def store_posterior_samples(self):
        """Samples from the Posterior objects and sets raw_counts
        """

        for m, post in self.posteriors.items():
            pp_counts, original = post.sequential().generate(
                n_samples=self.n_samples, batch_size=self.batch_size
            )
            self.posterior_predictive_samples[m] = pp_counts
        self.raw_counts = original

    def coeff_of_variation(self, cell_wise: bool = True):
        """Calculate the coefficient of variation

        Args:
            cell_wise (bool, optional): Calculate for each cell across genes. Defaults to True.
                                        If False, calculate for each gene across cells.
        """
        axis = 1 if cell_wise is True else 0
        identifier = "cv_cell" if cell_wise is True else "cv_gene"
        df = pd.DataFrame()
        for m, samples in self.posterior_predictive_samples.items():
            cv = np.nanmean(
                np.std(samples, axis=axis) / np.mean(samples, axis=axis), axis=-1
            )

            df[m] = cv.ravel()

        df["raw"] = np.std(self.raw_counts, axis=axis) / np.mean(
            self.raw_counts, axis=axis
        )
        df["raw"] = np.nan_to_num(df["raw"])

        self.metrics[identifier] = df

    def mean(self, cell_wise: bool = False):
        """Calculate the mean across cells in one gene (or vice versa). Reports the mean and std over samples per model

        Args:
            cell_wise (bool, optional): Calculate for each cell across genes. Defaults to True.
                                        If False, calculate for each gene across cells.
        """
        axis = 1 if cell_wise is True else 0
        identifier = "mean_cell" if cell_wise is True else "mean_gene"
        df = pd.DataFrame()
        for m, samples in self.posterior_predictive_samples.items():
            item = np.mean(samples, axis=axis)
            item_mean = np.nanmean(item, axis=-1)
            item_std = np.nanstd(item, axis=-1)
            # make all zeros have 0 cv
            df[m + "_mean"] = item_mean.ravel()
            df[m + "_std"] = item_std.ravel()

        df["raw"] = np.mean(self.raw_counts, axis=axis)
        df["raw"] = np.nan_to_num(df["raw"])

        self.metrics[identifier] = df

    def variance(self, cell_wise: bool = False):
        """Calculate the mean across cells in one gene (or vice versa). Reports the mean and std over samples per model

        Args:
            cell_wise (bool, optional): Calculate for each cell across genes. Defaults to True.
                                        If False, calculate for each gene across cells.
        """
        axis = 1 if cell_wise is True else 0
        identifier = "var_cell" if cell_wise is True else "var_gene"
        df = pd.DataFrame()
        for m, samples in self.posterior_predictive_samples.items():
            item = np.var(samples, axis=axis)
            item_mean = np.nanmean(item, axis=-1)
            item_std = np.nanstd(item, axis=-1)
            # make all zeros have 0 cv
            df[m + "_mean"] = item_mean.ravel()
            df[m + "_std"] = item_std.ravel()

        df["raw"] = np.var(self.raw_counts, axis=axis)
        df["raw"] = np.nan_to_num(df["raw"])

        self.metrics[identifier] = df

    def median_absolute_error(self, point_estimate="mean"):
        df = pd.DataFrame()
        for m, samples in self.posterior_predictive_samples.items():
            if point_estimate == "mean":
                point_sample = np.mean(samples, axis=-1)
            else:
                point_sample = np.median(samples, axis=-1)
            mad_gene = np.median(
                np.abs(
                    point_sample[:, : self.dataset.nb_genes]
                    - self.raw_counts[:, : self.dataset.nb_genes]
                )
            )
            # For samples with protein data
            if point_sample.shape[1] != self.dataset.nb_genes:
                mad_pro = np.median(
                    np.abs(
                        point_sample[:, self.dataset.nb_genes :]
                        - self.raw_counts[:, self.dataset.nb_genes :]
                    )
                )
            else:
                mad_pro = np.nan
            df[m] = [mad_gene, mad_pro]

        df.index = ["genes", "proteins"]
        self.metrics["mae"] = df

    def mean_squared_log_error(self, point_estimate="mean"):
        df = pd.DataFrame()
        for m, samples in self.posterior_predictive_samples.items():
            if point_estimate == "mean":
                point_sample = np.mean(samples, axis=-1)
            else:
                point_sample = np.mean(samples, axis=-1)
            mad_gene = np.mean(
                np.square(
                    np.log(point_sample[:, : self.dataset.nb_genes] + 1)
                    - np.log(self.raw_counts[:, : self.dataset.nb_genes] + 1)
                )
            )
            if point_sample.shape[1] != self.dataset.nb_genes:
                mad_pro = np.mean(
                    np.square(
                        np.log(point_sample[:, self.dataset.nb_genes :] + 1)
                        - np.log(self.raw_counts[:, self.dataset.nb_genes :] + 1)
                    )
                )
            else:
                mad_pro = np.nan
            df[m] = [mad_gene, mad_pro]

        df.index = ["genes", "proteins"]
        self.metrics["msle"] = df

    def dropout_ratio(self):
        """Fraction of zeros in raw_counts for a specific gene
        """
        df = pd.DataFrame()
        for m, samples in self.posterior_predictive_samples.items():
            dr = np.mean(np.mean(samples == 0, axis=0), axis=-1)
            df[m] = dr.ravel()

        df["raw"] = np.mean(np.mean(self.raw_counts == 0, axis=0), axis=-1)

        self.metrics["dropout_ratio"] = df

    def store_external_samples(self, samples: np.ndarray, key: str):
        """ Samples should be cells by genes by posterior predictive samples

        """

        self.posterior_predictive_samples[key] = samples

    def store_fa_samples(self, key="Factor Analysis", normalization="log", **kwargs):
        # reconstruction
        if normalization == "log":
            data = np.log(self.raw_counts + 1)
            key += " (Log)"
        elif normalization == "rate":
            lib_size_rna = self.raw_counts[:, : self.dataset.nb_genes].sum(axis=1)[
                :, np.newaxis
            ]

            data = np.log(
                10000 * self.raw_counts[:, : self.dataset.nb_genes] / lib_size_rna + 1
            )
            key += " (Rate)"
        elif normalization == "log_rate":
            lib_size_rna = self.raw_counts[:, : self.dataset.nb_genes].sum(axis=1)[
                :, np.newaxis
            ]

            data = np.log(
                10000 * self.raw_counts[:, : self.dataset.nb_genes] / lib_size_rna + 1
            )
            key += " (Log Rate)"
        else:
            data = self.raw_counts
        fa = FactorAnalysis(**kwargs)
        fa.fit(data)
        self.models[key] = fa

        # transform gives the posterior mean
        z_mean = fa.transform(data)
        Ih = np.eye(len(fa.components_))
        # W is n_components by n_features, code below from sklearn implementation
        Wpsi = fa.components_ / fa.noise_variance_
        z_cov = linalg.inv(Ih + np.dot(Wpsi, fa.components_.T))

        # sample z's
        z_samples = np.random.multivariate_normal(
            np.zeros(fa.n_components),
            cov=z_cov,
            size=(self.raw_counts.shape[0], self.n_samples),
        )
        # cells by n_components by posterior samples
        z_samples = np.swapaxes(z_samples, 1, 2)
        # add mean to all samples
        z_samples += z_mean[:, :, np.newaxis]

        x_samples = np.zeros(
            (self.raw_counts.shape[0], self.raw_counts.shape[1], self.n_samples)
        )
        for i in range(self.n_samples):
            x_mean = np.matmul(z_samples[:, :, i], fa.components_)
            x_sample = np.random.normal(x_mean, scale=np.sqrt(fa.noise_variance_))
            # add back feature means
            x_samples[:, :, i] = x_sample + fa.mean_

        reconstruction = x_samples

        if normalization == "log":
            reconstruction = np.exp(reconstruction - 1)
        if normalization == "rate":
            reconstruction = (
                lib_size_rna[:, :, np.newaxis]
                / 10000
                * reconstruction[:, : self.dataset.nb_genes],
            )
        if normalization == "log_rate":
            reconstruction = (
                lib_size_rna[:, :, np.newaxis]
                / 10000
                * np.exp(reconstruction[:, : self.dataset.nb_genes] - 1),
            )

        self.posterior_predictive_samples[key] = reconstruction

    def store_pca_samples(self, key="PCA", normalization="log", **kwargs):
        # reconstruction
        if normalization == "log":
            data = np.log(self.raw_counts + 1)
            key += " (Log)"
        elif normalization == "rate":
            lib_size_rna = self.raw_counts[:, : self.dataset.nb_genes].sum(axis=1)[
                :, np.newaxis
            ]

            data = np.log(
                10000 * self.raw_counts[:, : self.dataset.nb_genes] / lib_size_rna + 1
            )
            key += " (Rate)"
        elif normalization == "log_rate":
            lib_size_rna = self.raw_counts[:, : self.dataset.nb_genes].sum(axis=1)[
                :, np.newaxis
            ]

            data = np.log(
                10000 * self.raw_counts[:, : self.dataset.nb_genes] / lib_size_rna + 1
            )
            key += " (Log Rate)"
        else:
            data = self.raw_counts
        pca = PCA(**kwargs)
        pca.fit(data)
        self.models[key] = pca

        # Using Bishop notation section 12.2, M is comp x comp
        # W is fit using MLE, samples generated using posterior predictive
        M = (
            np.matmul(pca.components_, pca.components_.T)
            + np.identity(pca.n_components) * pca.noise_variance_
        )
        z_mean = np.matmul(
            np.matmul(linalg.inv(M), pca.components_), (self.raw_counts - pca.mean_).T
        ).T
        z_cov = linalg.inv(M) * pca.noise_variance_

        # sample z's
        z_samples = np.random.multivariate_normal(
            np.zeros(pca.n_components),
            cov=z_cov,
            size=(self.raw_counts.shape[0], self.n_samples),
        )
        # cells by n_components by posterior samples
        z_samples = np.swapaxes(z_samples, 1, 2)
        # add mean to all samples
        z_samples += z_mean[:, :, np.newaxis]

        x_samples = np.zeros(
            (self.raw_counts.shape[0], self.raw_counts.shape[1], self.n_samples)
        )
        for i in range(self.n_samples):
            x_mean = np.matmul(z_samples[:, :, i], pca.components_)
            x_sample = np.random.normal(x_mean, scale=np.sqrt(pca.noise_variance_))
            # add back feature means
            x_samples[:, :, i] = x_sample + pca.mean_

        reconstruction = x_samples

        if normalization == "log":
            reconstruction = np.clip(reconstruction, -1000, 20)
            reconstruction = np.exp(reconstruction - 1)
        if normalization == "rate":
            reconstruction = (
                lib_size_rna[:, :, np.newaxis]
                / 10000
                * reconstruction[:, : self.dataset.nb_genes],
            )
        if normalization == "log_rate":
            reconstruction = (
                lib_size_rna[:, :, np.newaxis]
                / 10000
                * np.exp(reconstruction[:, : self.dataset.nb_genes] - 1),
            )

        self.posterior_predictive_samples[key] = reconstruction

    def gene_gene_correlation(self, gene_indices=None, n_genes=1000):
        if gene_indices is not None:
            self.gene_set = np.random.choice(
                self.dataset.nb_genes, size=n_genes, replace=False
            )
        else:
            self.gene_set = gene_indices
        model_corrs = {}
        for m, samples in tqdm(self.posterior_predictive_samples.items()):
            correlation_matrix = np.zeros((n_genes, len(self.dataset.protein_names)))
            for i in range(self.n_samples):
                sample = StandardScaler().fit_transform(samples[:, :, i])
                gene_sample = sample[:, self.gene_set]
                correlation_matrix += np.matmul(gene_sample.T, gene_sample)
            correlation_matrix /= self.n_samples
            correlation_matrix /= self.raw_counts.shape[0] - 1
            model_corrs[m] = correlation_matrix.ravel()

        scaled_raw_counts = StandardScaler().fit_transform(self.raw_counts)
        scaled_genes = scaled_raw_counts[:, self.gene_set]
        raw_count_corr = np.matmul(scaled_genes.T, scaled_genes)
        raw_count_corr /= self.raw_counts.shape[0] - 1
        model_corrs["raw"] = raw_count_corr.ravel()

        model_corrs["gene_names1"] = (
            list(self.dataset.gene_names[self.gene_set]) * n_genes
        )
        model_corrs["gene_names2"] = np.repeat(
            self.dataset.gene_names[self.gene_set],
            len(self.dataset.gene_names[self.gene_set]),
        )

        df = pd.DataFrame.from_dict(model_corrs)
        self.metrics["all gene-gene correlations"] = df

    def calibration_error(self, confidence_intervals=None):
        if confidence_intervals is None:
            ps = [2.5, 5, 7.5, 10, 12.5, 15, 17.5, 82.5, 85, 87.5, 90, 92.5, 95, 97.5]
        else:
            ps = confidence_intervals
        reverse_ps = ps[::-1]
        model_cal = {}
        for m, samples in self.posterior_predictive_samples.items():
            percentiles = np.percentile(samples, ps, axis=2)
            reverse_percentiles = percentiles[::-1]
            cal_error_genes = 0
            for i, j, truth, reverse_truth in zip(
                percentiles, reverse_percentiles, ps, reverse_ps
            ):
                if truth > reverse_truth:
                    break
                true_width = (100 - truth * 2) / 100
                # For gene only model
                ci = np.logical_and(
                    self.raw_counts[:, : self.dataset.nb_genes] >= i,
                    self.raw_counts[:, : self.dataset.nb_genes] <= j,
                )
                pci_genes = np.mean(ci[:, : self.dataset.nb_genes])
                cal_error_genes += (pci_genes - true_width) ** 2
            model_cal[m] = {"genes": cal_error_genes}
        self.metrics["calibration"] = pd.DataFrame.from_dict(model_cal)
