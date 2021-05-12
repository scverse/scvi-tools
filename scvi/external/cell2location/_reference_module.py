import numpy as np
import pandas as pd
import pyro
import pyro.distributions as dist
import torch
from anndata import AnnData
from pyro.distributions import constraints
from pyro.distributions.transforms import SoftplusTransform
from pyro.infer.autoguide import init_to_mean
from pyro.nn import PyroModule
from scipy.sparse import csr_matrix, issparse
from torch.distributions import biject_to, transform_to

import scvi
from scvi import _CONSTANTS
from scvi.data._anndata import get_from_registry
from scvi.distributions._negative_binomial import _convert_mean_disp_to_counts_logits
from scvi.model.base import BaseModelClass
from scvi.module.base import PyroBaseModuleClass
from scvi.nn import one_hot

from ._base import AutoGuideMixinModule, PltExportMixin, TrainSampleMixin


def compute_cluster_averages(adata, labels, use_raw=True, layer=None):
    """
    Compute average expression of each gene in each cluster

    Parameters
    ----------
    adata
        AnnData object of reference single-cell dataset
    labels
        Name of adata.obs column containing cluster labels
    use_raw
        Use raw slow in adata?
    layer
        use layer in adata? provide layer name

    Returns
    -------
    pd.DataFrame of cluster average expression of each gene

    """

    if layer is not None:
        x = adata.layers["layer"]
        var_names = adata.var_names
    else:
        if not use_raw:
            x = adata.X
            var_names = adata.var_names
        else:
            if not adata.raw:
                raise ValueError(
                    "AnnData object has no raw data, change `use_raw=True, layer=None` or fix your object"
                )
            x = adata.raw.X
            var_names = adata.raw.var_names

    if sum(adata.obs.columns == labels) != 1:
        raise ValueError("cluster_col is absent in adata_ref.obs or not unique")

    all_clusters = np.unique(adata.obs[labels])
    averages_mat = np.zeros((1, x.shape[1]))

    for c in all_clusters:
        sparse_subset = csr_matrix(x[np.isin(adata.obs[labels], c), :])
        aver = sparse_subset.mean(0)
        averages_mat = np.concatenate((averages_mat, aver))
    averages_mat = averages_mat[1:, :].T
    averages_df = pd.DataFrame(data=averages_mat, index=var_names, columns=all_clusters)

    return averages_df


@biject_to.register(constraints.positive)
@transform_to.register(constraints.positive)
def _transform_to_positive(constraint):
    return SoftplusTransform()


class RegressionBackgroundDetectionTechModule(PyroModule):
    """
    Given cell type annotation for each cell, the corresponding reference cell type signatures :math:`g_{f,g}`,
    which represent the average mRNA count of each gene `g` in each cell type `f={1, .., F}`,
    are estimated from sc/snRNA-seq data using Negative Binomial regression,
    which allows to robustly combine data across technologies and batches.

    This model combines batches, and treats data :math:`D` as Negative Binomial distributed,
    given mean :math:`\mu` and overdispersion :math:`\alpha`:

    .. math::
        D_{c,g} \sim \mathtt{NB}(alpha=\alpha_{g}, mu=\mu_{c,g})
    .. math::
        \mu_{c,g} = (\mu_{f,g} + s_{e,g}) * y_c * y_{t,g}

    Which is equivalent to:

    .. math::
        D_{c,g} \sim \mathtt{Poisson}(\mathtt{Gamma}(\alpha_{f,g}, \alpha_{f,g} / \mu_{c,g}))

    Here, :math:`\mu_{f,g}` denotes average mRNA count in each cell type :math:`f` for each gene :math:`g`;
    :math:`y_c` denotes normalisation for each cell :math:`c` with a prior mean for each experiment :math:`e`,
        to account for RNA capture, sequencing depth.
    :math:`y_{t,g}` denotes per gene :math:`g` detection efficiency normalisation for each technology :math:`t`

    """

    def __init__(
        self,
        n_obs,
        n_vars,
        n_factors,
        n_batch,
        n_tech=None,
        batch_size=None,
        alpha_g_phi_hyp_prior={"alpha": 9.0, "beta": 3.0},
        gene_add_alpha_hyp_prior={"alpha": 9.0, "beta": 3.0},
        gene_add_mean_hyp_prior={
            "alpha": 1.0,
            "beta": 100.0,
        },
        detection_hyp_prior={"alpha": 200.0, "mean_alpha": 1.0, "mean_beta": 1.0},
        gene_tech_prior={"mean": 1, "alpha": 200},
        init_vals=None,
    ):
        """

        Parameters
        ----------
        n_obs
        n_vars
        n_factors
        n_batch
        n_tech
        batch_size
        alpha_g_phi_hyp_prior
        gene_add_alpha_hyp_prior
        gene_add_mean_hyp_prior
        detection_hyp_prior
        gene_tech_prior
        use_average_as_initial_value
        """

        ############# Initialise parameters ################
        super().__init__()

        self.n_obs = n_obs
        self.n_vars = n_vars
        self.n_factors = n_factors
        self.n_batch = n_batch
        self.n_tech = n_tech
        self.batch_size = batch_size

        self.alpha_g_phi_hyp_prior = alpha_g_phi_hyp_prior
        self.gene_add_alpha_hyp_prior = gene_add_alpha_hyp_prior
        self.gene_add_mean_hyp_prior = gene_add_mean_hyp_prior
        self.detection_hyp_prior = detection_hyp_prior
        self.gene_tech_prior = gene_tech_prior

        if (init_vals is not None) & (type(init_vals) is dict):
            self.init_vals = init_vals

        self.register_buffer(
            "detection_hyp_prior_alpha",
            torch.tensor(self.detection_hyp_prior["alpha"]),
        )
        self.register_buffer(
            "detection_mean_hyp_prior_alpha",
            torch.tensor(self.detection_hyp_prior["mean_alpha"]),
        )
        self.register_buffer(
            "detection_mean_hyp_prior_beta",
            torch.tensor(self.detection_hyp_prior["mean_beta"]),
        )
        self.register_buffer(
            "gene_tech_prior_alpha",
            torch.tensor(self.gene_tech_prior["alpha"]),
        )
        self.register_buffer(
            "gene_tech_prior_beta",
            torch.tensor(self.gene_tech_prior["alpha"] / self.gene_tech_prior["mean"]),
        )

        self.register_buffer(
            "alpha_g_phi_hyp_prior_alpha",
            torch.tensor(self.alpha_g_phi_hyp_prior["alpha"]),
        )
        self.register_buffer(
            "alpha_g_phi_hyp_prior_beta",
            torch.tensor(self.alpha_g_phi_hyp_prior["beta"]),
        )
        self.register_buffer(
            "gene_add_alpha_hyp_prior_alpha",
            torch.tensor(self.gene_add_alpha_hyp_prior["alpha"]),
        )
        self.register_buffer(
            "gene_add_alpha_hyp_prior_beta",
            torch.tensor(self.gene_add_alpha_hyp_prior["beta"]),
        )
        self.register_buffer(
            "gene_add_mean_hyp_prior_alpha",
            torch.tensor(self.gene_add_mean_hyp_prior["alpha"]),
        )
        self.register_buffer(
            "gene_add_mean_hyp_prior_beta",
            torch.tensor(self.gene_add_mean_hyp_prior["beta"]),
        )

        self.register_buffer("ones", torch.ones((1, 1)))
        self.register_buffer("eps", torch.tensor(1e-8))

    ############# Define the model ################
    @staticmethod
    def _get_fn_args_from_batch(tensor_dict):
        x_data = tensor_dict[_CONSTANTS.X_KEY]
        ind_x = tensor_dict["ind_x"].long().squeeze()
        batch_index = tensor_dict[_CONSTANTS.BATCH_KEY]
        label_index = tensor_dict[_CONSTANTS.LABELS_KEY]
        return (x_data, ind_x, batch_index, label_index), {}

    def _get_fn_args_full_data(self, adata):
        x_data = get_from_registry(adata, _CONSTANTS.X_KEY)
        if issparse(x_data):
            x_data = np.array(x_data.toarray())
        x_data = torch.tensor(x_data.astype("float32"))
        ind_x = torch.tensor(get_from_registry(adata, "ind_x"))
        batch_index = torch.tensor(get_from_registry(adata, _CONSTANTS.BATCH_KEY))
        label_index = torch.tensor(get_from_registry(adata, _CONSTANTS.LABELS_KEY))
        return (x_data, ind_x, batch_index, label_index), {}

    def create_plates(self, x_data, idx, batch_index, label_index):
        return pyro.plate("obs_plate", size=self.n_obs, dim=-2, subsample=idx)

    def list_obs_plate_vars(self):
        """Create a dictionary with the name of observation/minibatch plate,
        indexes of model args to provide to encoder,
        variable names that belong to the observation plate
        and the number of dimensions in non-plate axis of each variable"""

        return {
            "name": "obs_plate",
            "in": [0],  # expression data + (optional) batch index
            "sites": {
                "detection_y_c": 1,
            },
        }

    def forward(self, x_data, idx, batch_index, label_index):

        obs2sample = one_hot(batch_index, self.n_batch)
        obs2label = one_hot(label_index, self.n_factors)
        obs_plate = self.create_plates(x_data, idx, batch_index, label_index)

        # =====================Per-cluster average mRNA count ======================= #
        # \mu_{f,g}
        per_cluster_mu_fg = pyro.sample(
            "per_cluster_mu_fg",
            dist.Gamma(self.ones, self.ones)
            .expand([self.n_factors, self.n_vars])
            .to_event(2),
        )

        # =====================Gene-specific multiplicative component ======================= #
        # `y_{t, g}` per gene multiplicative effect that explains the difference in sensitivity betweeen genes
        if self.n_tech is not None:
            detection_tech_gene_tg = pyro.sample(
                "detection_tech_gene_tg",
                dist.Gamma(
                    self.ones * self.gene_tech_prior_alpha,
                    self.ones * self.gene_tech_prior_beta,
                )
                .expand([self.n_tech, self.n_vars])
                .to_event(2),
            )

        # =====================Cell-specific detection efficiency ======================= #
        # y_c with hierarchical mean prior
        detection_mean_y_e = pyro.sample(
            "detection_mean_y_e",
            dist.Gamma(
                self.ones * self.detection_mean_hyp_prior_alpha,
                self.ones * self.detection_mean_hyp_prior_beta,
            )
            .expand([self.n_batch, 1])
            .to_event(2),
        )

        beta = (self.ones * self.detection_hyp_prior_alpha) / (
            obs2sample @ detection_mean_y_e
        )
        with obs_plate:
            detection_y_c = pyro.sample(
                "detection_y_c",
                dist.Gamma(self.ones * self.detection_hyp_prior_alpha, beta),
            )  # (self.n_obs, 1)

        # =====================Gene-specific additive component ======================= #
        # s_{e,g} accounting for background, free-floating RNA
        s_g_gene_add_alpha_hyp = pyro.sample(
            "s_g_gene_add_alpha_hyp",
            dist.Gamma(
                self.gene_add_alpha_hyp_prior_alpha, self.gene_add_alpha_hyp_prior_beta
            ),
        )
        s_g_gene_add_mean = pyro.sample(
            "s_g_gene_add_mean",
            dist.Gamma(
                self.gene_add_mean_hyp_prior_alpha,
                self.gene_add_mean_hyp_prior_beta,
            )
            .expand([self.n_batch, 1])
            .to_event(2),
        )  # (self.n_batch)
        s_g_gene_add_alpha_e_inv = pyro.sample(
            "s_g_gene_add_alpha_e_inv",
            dist.Exponential(s_g_gene_add_alpha_hyp)
            .expand([self.n_batch, 1])
            .to_event(2),
        )  # (self.n_batch)
        s_g_gene_add_alpha_e = self.ones / s_g_gene_add_alpha_e_inv.pow(2)

        s_g_gene_add = pyro.sample(
            "s_g_gene_add",
            dist.Gamma(s_g_gene_add_alpha_e, s_g_gene_add_alpha_e / s_g_gene_add_mean)
            .expand([self.n_batch, self.n_vars])
            .to_event(2),
        )  # (self.n_batch, n_vars)

        # =====================Gene-specific overdispersion ======================= #
        alpha_g_phi_hyp = pyro.sample(
            "alpha_g_phi_hyp",
            dist.Gamma(
                self.alpha_g_phi_hyp_prior_alpha, self.alpha_g_phi_hyp_prior_beta
            ),
        )
        alpha_g_inverse = pyro.sample(
            "alpha_g_inverse",
            dist.Exponential(alpha_g_phi_hyp).expand([1, self.n_vars]).to_event(2),
        )  # (self.n_batch or 1, self.n_vars)

        # =====================Expected expression ======================= #

        # overdispersion
        alpha = self.ones / alpha_g_inverse.pow(2)
        # biological expression
        mu_biol = (
            obs2label @ per_cluster_mu_fg  # contaminating RNA
            + obs2sample @ s_g_gene_add
        ) * detection_y_c  # cell-specific normalisation
        if self.n_tech is not None:
            mu_biol = mu_biol * detection_tech_gene_tg  # gene-specific normalisation
        total_count, logits = _convert_mean_disp_to_counts_logits(
            mu_biol, alpha, eps=self.eps
        )

        # =====================DATA likelihood ======================= #
        # Likelihood (sampling distribution) of data_target & add overdispersion via NegativeBinomial
        with obs_plate:
            pyro.sample(
                "data_target",
                dist.NegativeBinomial(total_count=total_count, logits=logits),
                # NegativeBinomial(mu=mu_biol, theta=alpha_biol),
                obs=x_data,
            )

    # =====================Other functions======================= #
    def compute_expected(self, samples, adata):
        r"""Compute expected expression of each gene in each cell. Useful for evaluating how well
        the model learned expression pattern of all genes in the data.
        """
        obs2sample = get_from_registry(adata, _CONSTANTS.BATCH_KEY)
        obs2sample = pd.get_dummies(obs2sample.flatten())
        obs2label = get_from_registry(adata, _CONSTANTS.LABELS_KEY)
        obs2label = pd.get_dummies(obs2label.flatten())

        alpha = 1 / np.power(samples["alpha_g_inverse"], 2)

        mu = (
            np.dot(obs2label, samples["per_cluster_mu_fg"])
            + np.dot(obs2sample, samples["s_g_gene_add"])
        ) * samples["detection_y_c"]
        if self.n_tech is not None:
            mu = mu * samples["detection_tech_gene_tg"]

        return {"mu": mu, "alpha": alpha}

    def compute_expected_subset(self, samples, adata, fact_ind, cell_ind):
        r"""Compute expected expression of each gene in each cell that comes from
        a subset of factors (cell types) or cells.
        Useful for evaluating how well the model learned expression pattern of all genes in the data.
        """
        obs2sample = get_from_registry(adata, _CONSTANTS.BATCH_KEY)
        obs2sample = pd.get_dummies(obs2sample.flatten())
        obs2label = get_from_registry(adata, _CONSTANTS.LABELS_KEY)
        obs2label = pd.get_dummies(obs2label.flatten())

        alpha = 1 / np.power(samples["alpha_g_inverse"], 2)

        mu = (
            np.dot(
                obs2label[cell_ind, fact_ind], samples["per_cluster_mu_fg"][fact_ind, :]
            )
            + np.dot(obs2sample[cell_ind, :], samples["s_g_gene_add"])
        ) * samples["detection_y_c"]
        if self.n_tech is not None:
            mu = mu * samples["detection_tech_gene_tg"]

        return {"mu": mu, "alpha": alpha}

    def normalise_by_sample_scaling(self):
        r"""Normalise expression data by inferred technical variables (compute pearson residuals)."""
        None


class RegressionBaseModule(PyroBaseModuleClass, AutoGuideMixinModule):
    def __init__(
        self,
        model,
        amortised: bool = False,
        single_encoder: bool = True,
        encoder_kwargs=None,
        data_transform="log1p",
        **kwargs,
    ):
        """
        Module class which defines AutoGuide given model. Supports multiple model architectures.

        Parameters
        ----------
        amortised
            boolean, use a Neural Network to approximate posterior distribution of location-specific (local) parameters?
        encoder_kwargs
            arguments for Neural Network construction (scvi.nn.FCLayers)
        kwargs
            arguments for specific model class - e.g. number of genes, values of the prior distribution
        """
        super().__init__()
        self.hist = []

        self._model = model(**kwargs)
        self._amortised = amortised

        self._guide = self._create_autoguide(
            model=self.model,
            amortised=self.is_amortised,
            encoder_kwargs=encoder_kwargs,
            data_transform=data_transform,
            single_encoder=single_encoder,
            init_loc_fn=init_to_mean,
        )

        self._get_fn_args_from_batch = self._model._get_fn_args_from_batch

    @property
    def model(self):
        return self._model

    @property
    def guide(self):
        return self._guide

    @property
    def is_amortised(self):
        return self._amortised


class RegressionModel(TrainSampleMixin, BaseModelClass, PltExportMixin):
    """
    Model which estimates per cluster average mRNA count account for batch effects. User-end model class.

    https://github.com/BayraktarLab/cell2location

    Parameters
    ----------
    adata
        single-cell AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    use_gpu
        Use the GPU?
    **model_kwargs
        Keyword args for :class:`~scvi.external.LocationModelLinearDependentWMultiExperimentModel`

    Examples
    --------
    TODO add example
    >>>
    """

    def __init__(
        self,
        adata: AnnData,
        batch_size=None,
        model_class=None,
        **model_kwargs,
    ):
        # add index for each cell (provided to pyro plate for correct minibatching)
        adata.obs["_indices"] = np.arange(adata.n_obs).astype("int64")
        scvi.data.register_tensor_from_anndata(
            adata,
            registry_key="ind_x",
            adata_attr_name="obs",
            adata_key_name="_indices",
        )

        super().__init__(adata)

        if model_class is None:
            model_class = RegressionBackgroundDetectionTechModule

        self.batch_size = batch_size
        self.n_factors_ = self.summary_stats["n_labels"]
        self.factor_names_ = self.adata.uns["_scvi"]["categorical_mappings"][
            "_scvi_labels"
        ]["mapping"]

        self.module = RegressionBaseModule(
            model=model_class,
            n_obs=self.summary_stats["n_cells"],
            n_vars=self.summary_stats["n_vars"],
            n_factors=self.n_factors_,
            n_batch=self.summary_stats["n_batch"],
            **model_kwargs,
        )
        self._model_summary_string = f'RegressionBackgroundDetectionTech model with the following params: \nn_factors: {self.n_factors_} \nn_batch: {self.summary_stats["n_batch"]} '
        self.init_params_ = self._get_init_params(locals())
