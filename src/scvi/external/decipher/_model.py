from collections.abc import Sequence

import numpy as np
import pyro
import torch
import torch.nn.functional as F
from anndata import AnnData

from scvi._constants import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import LayerField
from scvi.model.base import BaseModelClass, PyroSviTrainMixin
from scvi.train import PyroTrainingPlan
from scvi.utils import setup_anndata_dsp

from ._module import DecipherPyroModule
from ._trainingplan import DecipherTrainingPlan
from .utils._trajectory import Trajectory


class Decipher(PyroSviTrainMixin, BaseModelClass):
    """Decipher model for single-cell data analysis :cite:p:`Nazaret24`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via
        :meth:`~scvi.model.Decipher.setup_anndata`.
    dim_v
        Dimension of the interpretable latent space v.
    dim_z
        Dimension of the intermediate latent space z.
    layers_v_to_z
        Hidden layer sizes for the v to z decoder network.
    layers_z_to_x
        Hidden layer sizes for the z to x decoder network.
    beta
        Regularization parameter for the KL divergence.
    """

    _module_cls = DecipherPyroModule
    _training_plan_cls = DecipherTrainingPlan

    def __init__(self, adata: AnnData, **kwargs):
        pyro.clear_param_store()

        super().__init__(adata)

        dim_genes = self.summary_stats.n_vars

        self.module = self._module_cls(
            dim_genes,
            **kwargs,
        )

        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        **kwargs,
    ) -> AnnData | None:
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    def train(
        self,
        max_epochs: int | None = None,
        accelerator: str = "auto",
        device: int | str = "auto",
        train_size: float = 0.9,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        early_stopping: bool = False,
        training_plan: PyroTrainingPlan | None = None,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **trainer_kwargs,
    ):
        if "early_stopping_monitor" not in trainer_kwargs:
            trainer_kwargs["early_stopping_monitor"] = "nll_validation"
        datasplitter_kwargs = datasplitter_kwargs or {}
        if "drop_last" not in datasplitter_kwargs:
            datasplitter_kwargs["drop_last"] = True
        super().train(
            max_epochs=max_epochs,
            accelerator=accelerator,
            device=device,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            batch_size=batch_size,
            early_stopping=early_stopping,
            plan_kwargs=plan_kwargs,
            training_plan=training_plan,
            datasplitter_kwargs=datasplitter_kwargs,
            **trainer_kwargs,
        )

    def get_latent_representation(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        give_z: bool = False,
    ) -> np.ndarray:
        """Get the latent representation of the data.

        Parameters
        ----------
        adata
            AnnData object with the data to get the latent representation of.
        indices
            Indices of the data to get the latent representation of.
        batch_size
            Batch size to use for the data loader.
        give_z
            Whether to return the intermediate latent space z or the top-level
            latent space v.
        """
        self._check_if_trained(warn=False)
        adata = self._validate_anndata(adata)

        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)
        latent_locs = []
        for tensors in scdl:
            x = tensors[REGISTRY_KEYS.X_KEY]
            x = torch.log1p(x)
            x = x.to(self.module.device)
            z_loc, _ = self.module.encoder_x_to_z(x)
            if give_z:
                latent_locs.append(z_loc)
            else:
                v_loc, _ = self.module.encoder_zx_to_v(torch.cat([z_loc, x], dim=-1))
                latent_locs.append(v_loc)
        return torch.cat(latent_locs).detach().cpu().numpy()

    def compute_imputed_gene_expression(
        self,
        adata: AnnData | None = None,
        indices: Sequence[int] | None = None,
        batch_size: int | None = None,
        compute_covariances: bool = False,
        v_obsm_key: str | None = None,
        z_obsm_key: str | None = None,
    ) -> np.ndarray | tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Impute gene expression from the decipher model.

        Parameters
        ----------
        adata
            The annotated data matrix.
        indices
            Indices of the data to get the latent representation of.
        batch_size
            Batch size to use for the data loader.
        compute_covariances
            Whether to compute the covariances between the Decipher v and each gene.
        v_obsm_key
            Key in `adata.obsm` to use for the Decipher v. Required if
            `compute_covariances` is True.
        z_obsm_key
            Key in `adata.obsm` to use for the Decipher z. Required if
            `compute_covariances` is True.

        Returns
        -------
        The imputed gene expression, and the covariances between the Decipher v and each gene
        if `compute_covariances` is True.
        """
        if compute_covariances and (v_obsm_key is None or z_obsm_key is None):
            raise ValueError(
                "`v_obsm_key` and `z_obsm_key` must be provided if `compute_covariances` is True."
            )

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(adata=adata, indices=indices, batch_size=batch_size)

        imputed_gene_expression_batches = []
        for tensors in scdl:
            x = tensors[REGISTRY_KEYS.X_KEY]
            z_loc, _, _, _ = self.module.guide(x)
            mu = self.module.decoder_z_to_x(z_loc)
            mu = F.softmax(mu, dim=-1)
            library_size = x.sum(axis=-1, keepdim=True)
            imputed_gene_expr = (library_size * mu.detach().cpu()).numpy()
            imputed_gene_expression_batches.append(imputed_gene_expr)
        imputed_gene_expression = np.concatenate(imputed_gene_expression_batches, axis=0)

        if compute_covariances:
            G = imputed_gene_expression.shape[1]
            v_gene_covariance = np.cov(
                imputed_gene_expression,
                y=adata.obsm[v_obsm_key],
                rowvar=False,
            )[:G, G:]
            z_gene_covariance = np.cov(
                imputed_gene_expression,
                y=adata.obsm[z_obsm_key],
                rowvar=False,
            )[:G, G:]
            return imputed_gene_expression, v_gene_covariance, z_gene_covariance

        return imputed_gene_expression

    @staticmethod
    def compute_decipher_time(
        adata: AnnData,
        cluster_obs_key: str,
        trajectory: Trajectory,
        n_neighbors: int = 10,
    ) -> np.ndarray:
        """Compute the decipher time for each cell, based on the inferred trajectories.

        The decipher time is computed by KNN regression of the cells'
        decipher v on the trajectories.

        Parameters
        ----------
        adata : AnnData
            The annotated data matrix.
        cluster_obs_key : str
            The key in adata.obs containing cluster assignments.
        trajectory : Trajectory
            A :class:`~scvi.external.decipher.utils.Trajectory` object containing
            the trajectory information.
        n_neighbors : int
            The number of neighbors to use for the KNN regression.

        Returns
        -------
        The decipher time of each cell.
        """
        try:
            from sklearn.neighbors import KNeighborsRegressor
        except ImportError as err:
            raise ImportError("Please install scikit-learn -- `pip install scikit-learn`") from err

        decipher_time = np.full(adata.n_obs, np.nan)

        knn = KNeighborsRegressor(n_neighbors=n_neighbors)
        knn.fit(trajectory.trajectory_latent, trajectory.trajectory_time)
        is_on_trajectory = adata.obs[cluster_obs_key].isin(trajectory.cluster_ids)
        cells_on_trajectory_idx = np.where(is_on_trajectory)[0]

        decipher_time[cells_on_trajectory_idx] = knn.predict(
            adata.obsm[trajectory.rep_key][cells_on_trajectory_idx]
        )
        return decipher_time

    def compute_gene_patterns(
        self,
        adata: AnnData,
        trajectory: Trajectory,
        l_scale: float = 10_000,
        n_samples: int = 100,
    ) -> dict[str, np.ndarray]:
        """Compute the gene patterns for a trajectory.

        The trajectory's points are sent through the decoders, thus defining distributions over the
        gene expression. The gene patterns are computed by sampling from these distribution.

        Parameters
        ----------
        adata : AnnData
            The annotated data matrix.
        trajectory : Trajectory
            A :class:`~scvi.external.decipher.utils.Trajectory` object containing
            the trajectory information.
        l_scale : float
            The library size scaling factor.
        n_samples : int
            The number of samples to draw from the decoder to compute the gene pattern statistics.

        Returns
        -------
        The gene patterns for the trajectory.
        Dictionary keys:
            - `mean`: the mean gene expression pattern
            - `q25`: the 25% quantile of the gene expression pattern
            - `q75`: the 75% quantile of the gene expression pattern
            - `times`: the times of the trajectory
        """
        adata = self._validate_anndata(adata)

        t_points = trajectory.trajectory_latent
        t_times = trajectory.trajectory_time

        t_points = torch.FloatTensor(t_points).to(self.module.device)
        z_mean, z_scale = self.module.decoder_v_to_z(t_points)
        z_scale = F.softplus(z_scale)

        z_samples = torch.distributions.Normal(z_mean, z_scale).sample(sample_shape=(n_samples,))

        gene_patterns = {}
        gene_patterns["mean"] = (
            F.softmax(self.module.decoder_z_to_x(z_mean), dim=-1).detach().cpu().numpy() * l_scale
        )

        gene_expression_samples = (
            F.softmax(self.module.decoder_z_to_x(z_samples), dim=-1).detach().cpu().numpy()
            * l_scale
        )
        gene_patterns["q25"] = np.quantile(gene_expression_samples, 0.25, axis=0)
        gene_patterns["q75"] = np.quantile(gene_expression_samples, 0.75, axis=0)

        gene_patterns["times"] = t_times

        return gene_patterns
