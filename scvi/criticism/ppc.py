import numpy as np
import pandas as pd
from sklearn.decomposition import FactorAnalysis
from typing import Dict, Iterable, List, Tuple, Union, Optional
from scvi.inference import Posterior
import logging

logger = logging.getLogger(__name__)


class PosteriorPredictiveCheck:
    """Generic class for posterior predictive checks (ppcs) on fitted models.

    This class is scVI's bases PPC class. It gives access to standard metrics, as
    well as building reconstructions from simpler models like Factor Analysis.
    """

    def __init__(self, posteriors_dict: Dict[str, Posterior], n_samples: int = 10):
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

        self.store_posterior_samples()

    def store_posterior_samples(self):
        """Samples from the Posterior objects and sets raw_counts
        """

        for m, post in self.posteriors.items():
            pp_counts, original = post.sequential().generate(n_samples=self.n_samples)
            self.posterior_predictive_samples[m] = pp_counts
        self.raw_counts = original

    def coeff_of_variation(self, cell_wise: bool = True):
        """Calculate the coefficient of variation
        
        Args:
            cell_wise (bool, optional): Calculate over cells. Defaults to True. If False, calculate over genes
        """
        axis = 0 if cell_wise is True else 1
        identifier = "cv_cell" if cell_wise is True else "cv_gene"
        df = pd.DataFrame()
        for m, samples in self.posterior_predictive_samples.items():
            cv = np.mean(
                np.std(samples, axis=axis) / np.mean(samples, axis=axis), axis=-1
            )
            df[m] = cv.ravel()

        df["raw"] = np.std(self.raw_counts, axis=axis) / np.mean(
            self.raw_counts, axis=axis
        )

        self.metrics[identifier] = df

    def dropout_ratio(self):
        """Fraction of zeros in raw_counts for a specific gene
        """
        df = pd.DataFrame()
        for m, samples in self.posterior_predictive_samples.items():
            dr = np.mean(np.mean(samples == 0, axis=0), axis=-1)
            df[m] = dr.ravel()

        df["raw"] = np.mean(np.mean(self.raw_counts == 0, axis=0), axis=-1)

        self.metrics["dropout_ratio"] = df

    def store_fa_samples(self, key="fa", **kwargs):
        self._fit_factor_analysis(key=key, **kwargs)
        # reconstruction
        fa = self.models["fa"]
        factors = fa.transform(self.raw_counts)
        comps = fa.components_
        cov = fa.get_covariance()
        recon_fa_mean = np.matmul(factors, comps)[:, :, np.newaxis]
        feature_means = np.mean(self.dataset.X, axis=0).A.ravel()[
            np.newaxis, :, np.newaxis
        ]
        samples = np.random.multivariate_normal(
            np.zeros(self.raw_counts.shape[1]),
            cov=cov,
            size=(self.dataset.X.shape[0], self.n_samples),
        )
        # cells by genes by posterior samples
        samples = np.swapaxes(samples, 1, 2)
        reconstruction = samples + recon_fa_mean + feature_means

        self.posterior_predictive_samples[key] = reconstruction

    def _fit_factor_analysis(self, key="fa", **kwargs):

        fa = FactorAnalysis(**kwargs)
        fa.fit(self.raw_counts)
        self.models[key] = fa
