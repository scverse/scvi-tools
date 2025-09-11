from __future__ import annotations

import logging
import warnings
from typing import TYPE_CHECKING

import numpy as np
import torch
import xarray as xr
from torch.distributions import Independent
from tqdm import tqdm

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager, fields
from scvi.data._utils import _validate_adata_dataloader_input
from scvi.external.mrvi._types import MRVIReduction
from scvi.external.mrvi_torch._module import TorchMRVAE
from scvi.external.mrvi_torch._utils import rowwise_max_excluding_diagonal
from scvi.model.base import BaseModelClass, UnsupervisedTrainingMixin, VAEMixin
from scvi.utils import setup_anndata_dsp

if TYPE_CHECKING:
    from collections.abc import Iterator
    from typing import Literal

    import numpy.typing as npt
    from anndata import AnnData
    from torch import Tensor
    from torch.distributions import Distribution

logger = logging.getLogger(__name__)


DEFAULT_TRAIN_KWARGS = {
    "max_epochs": 100,
    "early_stopping": True,
    "early_stopping_patience": 15,
    "check_val_every_n_epoch": 1,
    "batch_size": 256,
    "train_size": 0.9,
    "gradient_clip_val": 40,
    "plan_kwargs": {
        "lr": 2e-3,
        "n_epochs_kl_warmup": 20,
        "eps": 1e-8,
        "weight_decay": 1e-8,
    },
}


class TorchMRVI(VAEMixin, UnsupervisedTrainingMixin, BaseModelClass):
    """Multi-resolution Variational Inference (MrVI) :cite:p:`Boyeau24`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.external.TorchMRVI.setup_anndata`.
    n_latent
        Dimensionality of the latent space for ``z``.
    n_latent_u
        Dimensionality of the latent space for ``u``.
    encoder_n_hidden
        Number of nodes per hidden layer in the encoder.
    encoder_n_layers
        Number of hidden layers in the encoder.
    z_u_prior
        Whether to use a prior for ``z_u``.
    z_u_prior_scale
        Scale of the prior for the difference between ``z`` and ``u``.
    u_prior_scale
        Scale of the prior for ``u``.
    u_prior_mixture
        Whether to use a mixture model for the ``u`` prior.
    u_prior_mixture_k
        Number of components in the mixture model for the ``u`` prior.
    learn_z_u_prior_scale
        Whether to learn the scale of the ``z`` and ``u`` difference prior during training.
    laplace_scale
        Scale parameter for the Laplace distribution in the decoder.
    scale_observations
        Whether to scale loss by the number of observations per sample.
    px_kwargs
        Keyword args for :class:`~scvi.external.mrvi_torch._module.DecoderZXAttention`.
    qz_kwargs
        Keyword args for :class:`~scvi.external.mrvi_torch._module.EncoderUZ`.
    qu_kwargs
        Keyword args for :class:`~scvi.external.mrvi_torch._module.EncoderXU`.

    Notes
    -----
    This implementation of MRVI is in PyTorch.
    This will become the default version in v1.5 for MRVI.

    See further usage examples in the following tutorial:

    1. :doc:`/tutorials/notebooks/scrna/MrVI_tutorial_torch`

    See the user guide for this model:

    1. :doc:`/user_guide/models/mrvi_jax`

    See Also
    --------
    :class:`~scvi.external.mrvi_torch.TorchMRVAE`
    """

    _module_cls = TorchMRVAE

    def __init__(self, adata: AnnData | None = None, registry: dict | None = None, **model_kwargs):
        super().__init__(adata, registry)

        warnings.warn(
            "You are using the Torch Version of MrVI, starting v1.5, This will become the "
            "default version of MrVI instead of the Jax backend version.",
            DeprecationWarning,
            stacklevel=settings.warnings_stacklevel,
        )

        n_batch = self.summary_stats.n_batch
        n_labels = self.summary_stats.n_labels
        n_sample = self.summary_stats.n_sample

        if adata is not None:
            self.update_sample_info(adata)
            self.sample_key = self.adata_manager.get_state_registry(
                REGISTRY_KEYS.SAMPLE_KEY
            ).original_key
            self.sample_order = self.adata_manager.get_state_registry(
                REGISTRY_KEYS.SAMPLE_KEY
            ).categorical_mapping
            self.n_obs_per_sample = torch.Tensor(
                adata.obs._scvi_sample.value_counts().sort_index().values
            )
        else:
            self.sample_key = self.registry["field_registries"][f"{REGISTRY_KEYS.SAMPLE_KEY}"][
                "state_registry"
            ]["original_key"]
            self.sample_order = self.registry["field_registries"][f"{REGISTRY_KEYS.SAMPLE_KEY}"][
                "state_registry"
            ]["categorical_mapping"]
            self.n_obs_per_sample = self.registry["field_registries"][
                f"{REGISTRY_KEYS.SAMPLE_KEY}"
            ]["n_obs_per_sample"]["n_obs_per_sample"]
        self.backend = "torch"

        self.module = self._module_cls(
            n_input=self.summary_stats.n_vars,
            n_sample=n_sample,
            n_batch=n_batch,
            n_labels=n_labels,
            n_obs_per_sample=self.n_obs_per_sample,
            **model_kwargs,
        )
        self.init_params_ = self._get_init_params(locals())

    def to_device(self, device):
        pass

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        sample_key: str | None = None,
        batch_key: str | None = None,
        labels_key: str | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        %(param_sample_key)s
        %(param_batch_key)s
        %(param_labels_key)s
        **kwargs
            Additional keyword arguments passed into
            :meth:`~scvi.data.AnnDataManager.register_fields`.
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        # Add index for batched computation of local statistics.
        adata.obs["_indices"] = np.arange(adata.n_obs).astype(int)
        anndata_fields = [
            fields.LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            fields.CategoricalObsField(REGISTRY_KEYS.SAMPLE_KEY, sample_key),
            fields.CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            fields.CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
            fields.NumericalObsField(REGISTRY_KEYS.INDICES_KEY, "_indices"),
        ]

        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    def train(
        self,
        max_epochs: int | None = None,
        accelerator: str | None = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float | None = None,
        validation_size: float | None = None,
        batch_size: int = 128,
        early_stopping: bool = False,
        plan_kwargs: dict | None = None,
        **trainer_kwargs,
    ):
        """Train the model.

        Parameters
        ----------
        max_epochs
            Maximum number of epochs to train the model. The actual number of epochs may be less if
            early stopping is enabled. If ``None``, defaults to a heuristic based on
            :func:`~scvi.model.get_max_epochs_heuristic`.
        %(param_accelerator)s
        %(param_devices)s
        train_size
            Size of training set in the range ``[0.0, 1.0]``.
        validation_size
            Size of the validation set. If ``None``, defaults to ``1 - train_size``. If
            ``train_size + validation_size < 1``, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        early_stopping
            Perform early stopping. Additional arguments can be passed in through ``**kwargs``.
            See :class:`~scvi.train.Trainer` for further options.
        plan_kwargs
            Additional keyword arguments passed into :class:`~scvi.train.TrainingPlan`.
        **trainer_kwargs
            Additional keyword arguments passed into :class:`~scvi.train.Trainer`.
        """
        from copy import deepcopy

        train_kwargs = {
            "max_epochs": max_epochs,
            "accelerator": accelerator,
            "devices": devices,
            "train_size": train_size,
            "validation_size": validation_size,
            "batch_size": batch_size,
            "early_stopping": early_stopping,
            **trainer_kwargs,
        }
        train_kwargs = dict(deepcopy(DEFAULT_TRAIN_KWARGS), **train_kwargs)
        plan_kwargs = plan_kwargs or {}
        train_kwargs["plan_kwargs"] = dict(
            deepcopy(DEFAULT_TRAIN_KWARGS["plan_kwargs"]), **plan_kwargs
        )

        super().train(**train_kwargs)

    def get_latent_representation(
        self,
        adata: AnnData | None = None,
        indices: npt.ArrayLike | None = None,
        batch_size: int | None = None,
        use_mean: bool = True,
        give_z: bool = False,
        dataloader: Iterator[dict[str, Tensor | None]] | None = None,
    ) -> npt.NDArray:
        """Compute the latent representation of the data.

        Parameters
        ----------
        adata
            AnnData object to use. Defaults to the AnnData object used to initialize the model.
        indices
            Indices of cells to use.
        batch_size
            Batch size to use for computing the latent representation.
        use_mean
            Whether to use the mean of the distribution as the latent representation.
        give_z
            Whether to return the z latent representation or the u latent representation.
        dataloader
            An iterator over minibatches of data on which to compute the metric. The minibatches
            should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
            the model. If ``None``, a dataloader is created from ``adata``.

        Returns
        -------
        The latent representation of the data.
        """
        self._check_if_trained(warn=False)
        _validate_adata_dataloader_input(self, adata, dataloader)
        if dataloader is None:
            adata = self._validate_anndata(adata)
            dataloader = self._make_data_loader(
                adata=adata,
                indices=indices,
                batch_size=batch_size,
            )
        else:
            for param in [indices, batch_size]:
                if param is not None:
                    Warning(
                        f"Using {param} after custom Dataloader was initialize is redundant, "
                        f"please re-initialize with selected {param}",
                    )
        us = []
        zs = []

        for tensors in dataloader:
            inference_inputs = self.module._get_inference_input(tensors)
            inference_inputs["use_mean"] = use_mean
            outputs = self.module.inference(**inference_inputs)

            if give_z:
                zs.append(outputs["z"].detach().cpu())
            else:
                us.append(outputs["u"].detach().cpu())

        if give_z:
            return np.concatenate(zs, axis=0)
        else:
            return np.concatenate(us, axis=0)

    def compute_local_statistics(
        self,
        reductions: list[MRVIReduction],
        adata: AnnData | None = None,
        indices: npt.ArrayLike | None = None,
        batch_size: int | None = None,
        use_vmap: Literal["auto", True, False] = "auto",
        norm: str = "l2",
        mc_samples: int = 10,
    ) -> xr.Dataset:
        """Compute local statistics from counterfactual sample representations.

        Local statistics are reductions over either the local counterfactual latent representations
        or the resulting local sample distance matrices. For a large number of cells and/or
        samples, this method can avoid scalability issues by grouping over cell-level covariates.

        Parameters
        ----------
        reductions
            List of reductions to compute over local counterfactual sample representations.
        adata
            AnnData object to use.
        indices
            Indices of cells to use.
        batch_size
            Batch size to use for computing the local statistics.
        use_vmap
            Whether to use vmap to compute the local statistics. If "auto", vmap will be used if
            the number of samples is less than 500.
        norm
            Norm to use for computing the distances.
        mc_samples
            Number of Monte Carlo samples to use for computing the local statistics. Only applies
            if using sampled representations.
        """
        from functools import partial

        from scvi.external.mrvi._types import _parse_local_statistics_requirements

        use_vmap = use_vmap if use_vmap != "auto" else self.summary_stats.n_sample < 500

        if not reductions or len(reductions) == 0:
            raise ValueError("At least one reduction must be provided.")

        adata = self.adata if adata is None else adata
        self._check_if_trained(warn=False)
        # Hack to ensure new AnnDatas have indices.
        adata.obs["_indices"] = np.arange(adata.n_obs).astype(int)

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size, iter_ndarray=True
        )
        n_sample = self.summary_stats.n_sample

        reqs = _parse_local_statistics_requirements(reductions)

        def mapped_inference_fn(
            x: torch.Tensor,
            sample_index: torch.Tensor,
            cf_sample: torch.Tensor,
            use_mean: bool,
            mc_samples: int | None = None,
        ):
            def inference_fn(cf):
                return self.module.inference(
                    x, sample_index, cf_sample=cf, mc_samples=mc_samples, use_mean=use_mean
                )["z"]

            inference_partial = partial(
                self.module.inference, x, sample_index, mc_samples=mc_samples, use_mean=use_mean
            )

            if use_vmap:
                return torch.vmap(inference_fn, in_dims=0, out_dims=-2, randomness="different")(
                    cf_sample
                )
            else:

                def per_sample_inference_fn(cf_sample):
                    return inference_partial(cf_sample=cf_sample)["z"]

                return torch.stack(
                    [per_sample_inference_fn(cf_sample=cf_sample_i) for cf_sample_i in cf_sample]
                ).permute(1, 0, 2)  # (n_cells, n_samples, n_latent)

        ungrouped_data_arrs = {}
        grouped_data_arrs = {}
        for ur in reqs.ungrouped_reductions:
            ungrouped_data_arrs[ur.name] = []
        for gr in reqs.grouped_reductions:
            grouped_data_arrs[gr.name] = {}  # Will map group category to running group sum.

        for tensor in tqdm(scdl):
            indices = tensor[REGISTRY_KEYS.INDICES_KEY].astype(int).flatten()
            n_cells = tensor[REGISTRY_KEYS.X_KEY].shape[0]
            cf_sample = np.broadcast_to(np.arange(n_sample)[:, None, None], (n_sample, n_cells, 1))
            inf_inputs = self.module._get_inference_input(
                tensor,
            )
            if reqs.needs_mean_representations:
                try:
                    mean_zs_ = mapped_inference_fn(
                        x=torch.Tensor(inf_inputs["x"]),
                        sample_index=torch.Tensor(inf_inputs["sample_index"]),
                        cf_sample=torch.Tensor(cf_sample),
                        use_mean=True,
                    )
                except RuntimeError as e:
                    if use_vmap:
                        raise RuntimeError(
                            # TODO: update error message
                            "Out of memory. Try setting use_vmap=False."
                        ) from e
                    else:
                        raise e

                mean_zs = xr.DataArray(
                    mean_zs_.detach().cpu().numpy(),
                    dims=["cell_name", "sample", "latent_dim"],
                    coords={
                        "cell_name": adata.obs_names[indices].values,
                        "sample": self.sample_order,
                    },
                    name="sample_representations",
                )
            if reqs.needs_sampled_representations:
                sampled_zs_ = mapped_inference_fn(
                    x=torch.Tensor(inf_inputs["x"]),
                    sample_index=torch.Tensor(inf_inputs["sample_index"]),
                    cf_sample=torch.Tensor(cf_sample),
                    use_mean=False,
                    mc_samples=mc_samples,
                )  # (n_mc_samples, n_cells, n_samples, n_latent)
                sampled_zs_ = sampled_zs_.permute((1, 0, 2, 3))
                sampled_zs = xr.DataArray(
                    sampled_zs_.detach().cpu().numpy(),
                    dims=["cell_name", "mc_sample", "sample", "latent_dim"],
                    coords={
                        "cell_name": adata.obs_names[indices].values,
                        "sample": self.sample_order,
                    },
                    name="sample_representations",
                )

            if reqs.needs_mean_distances:
                mean_dists = self._compute_distances_from_representations(
                    mean_zs_, indices, norm=norm, return_numpy=True
                )

            if reqs.needs_sampled_distances or reqs.needs_normalized_distances:
                sampled_dists = self._compute_distances_from_representations(
                    sampled_zs_, indices, norm=norm, return_numpy=True
                )

                if reqs.needs_normalized_distances:
                    if norm != "l2":
                        raise ValueError(
                            f"Norm must be 'l2' when using normalized distances. Got {norm}."
                        )
                    (
                        normalization_means,
                        normalization_vars,
                    ) = self._compute_local_baseline_dists(
                        tensor, mc_samples=mc_samples
                    )  # both are shape (n_cells,)
                    normalization_means = normalization_means.reshape(-1, 1, 1, 1)
                    normalization_vars = normalization_vars.reshape(-1, 1, 1, 1)
                    normalized_dists = (
                        (sampled_dists - normalization_means) / (normalization_vars**0.5)
                    ).mean(dim="mc_sample")  # (n_cells, n_samples, n_samples)

            # Compute each reduction
            for r in reductions:
                if r.input == "mean_representations":
                    inputs = mean_zs
                elif r.input == "sampled_representations":
                    inputs = sampled_zs
                elif r.input == "mean_distances":
                    inputs = mean_dists
                elif r.input == "sampled_distances":
                    inputs = sampled_dists
                elif r.input == "normalized_distances":
                    inputs = normalized_dists
                else:
                    raise ValueError(f"Unknown reduction input: {r.input}")

                outputs = r.fn(inputs)

                if r.group_by is not None:
                    group_by = adata.obs[r.group_by].iloc[indices]
                    group_by_cats = group_by.unique()
                    for cat in group_by_cats:
                        cat_summed_outputs = outputs.sel(
                            cell_name=adata.obs_names[indices][group_by == cat].values
                        ).sum(dim="cell_name")
                        cat_summed_outputs = cat_summed_outputs.assign_coords(
                            {f"{r.group_by}_name": cat}
                        )
                        if cat not in grouped_data_arrs[r.name]:
                            grouped_data_arrs[r.name][cat] = cat_summed_outputs
                        else:
                            grouped_data_arrs[r.name][cat] += cat_summed_outputs
                else:
                    ungrouped_data_arrs[r.name].append(outputs)

        # Combine all outputs.
        final_data_arrs = {}
        for ur_name, ur_outputs in ungrouped_data_arrs.items():
            final_data_arrs[ur_name] = xr.concat(ur_outputs, dim="cell_name")

        for gr in reqs.grouped_reductions:
            group_by = adata.obs[gr.group_by]
            group_by_counts = group_by.value_counts()
            averaged_grouped_data_arrs = []
            for cat, count in group_by_counts.items():
                averaged_grouped_data_arrs.append(grouped_data_arrs[gr.name][cat] / count)
            final_data_arr = xr.concat(averaged_grouped_data_arrs, dim=f"{gr.group_by}_name")
            final_data_arrs[gr.name] = final_data_arr

        return xr.Dataset(data_vars=final_data_arrs)

    def _compute_local_baseline_dists(
        self, batch: dict, mc_samples: int = 250
    ) -> tuple[npt.NDArray, npt.NDArray]:
        """
        Approximate the distributions used as baselines for normalizing the local sample distances.

        Approximates the means and variances of the Euclidean distance between two samples of
        the z latent representation for the original sample for each cell in ``adata``.

        Reference: https://www.overleaf.com/read/mhdxcrknzxpm.

        Parameters
        ----------
        batch
            Batch of data to compute the local sample representation for.
        mc_samples
            Number of Monte Carlo samples to use for computing the local baseline distributions.
        """
        mc_samples_per_cell = (
            mc_samples * 2
        )  # need double for pairs of samples to compute distance between

        kwargs = self.module._get_inference_input(batch)
        kwargs.update({"use_mean": False, "mc_samples": mc_samples_per_cell})
        outputs = self.module.inference(**kwargs)

        z = outputs["z"]
        first_half_z, second_half_z = z[:mc_samples], z[mc_samples:]
        l2_dists = torch.sqrt(torch.sum((first_half_z - second_half_z) ** 2, axis=2)).T

        return (
            torch.mean(l2_dists, axis=1).detach().numpy(),
            torch.var(l2_dists, axis=1).detach().numpy(),
        )

    def _compute_distances_from_representations(
        self,
        reps: torch.Tensor,
        indices: torch.Tensor,
        norm: Literal["l2", "l1", "linf"] = "l2",
        return_numpy: bool = True,
    ) -> xr.DataArray:
        if norm not in ("l2", "l1", "linf"):
            raise ValueError(f"`norm` {norm} not supported")

        def _compute_distance(rep: torch.Tensor):
            delta_mat = torch.unsqueeze(rep, 0) - torch.unsqueeze(rep, 1)
            if norm == "l2":
                res = delta_mat**2
                res = torch.sqrt(res.sum(-1))
            elif norm == "l1":
                res = torch.abs(delta_mat).sum(-1)
            elif norm == "linf":
                res = torch.abs(delta_mat).max(-1)
            return res

        if reps.ndim == 3:
            dists = torch.vmap(_compute_distance)(reps)
            if return_numpy:
                dists = dists.detach().cpu().numpy()
            return xr.DataArray(
                dists,
                dims=["cell_name", "sample_x", "sample_y"],
                coords={
                    "cell_name": self.adata.obs_names[indices].values,
                    "sample_x": self.sample_order,
                    "sample_y": self.sample_order,
                },
                name="sample_distances",
            )
        else:
            # Case with sampled representations
            dists = torch.vmap(torch.vmap(_compute_distance))(reps)
            if return_numpy:
                dists = np.array(dists.detach().cpu().numpy())
            return xr.DataArray(
                dists,
                dims=["cell_name", "mc_sample", "sample_x", "sample_y"],
                coords={
                    "cell_name": self.adata.obs_names[indices].values,
                    "mc_sample": np.arange(reps.shape[1]),
                    "sample_x": self.sample_order,
                    "sample_y": self.sample_order,
                },
                name="sample_distances",
            )

    def get_local_sample_representation(
        self,
        adata: AnnData | None = None,
        indices: npt.ArrayLike | None = None,
        batch_size: int = 256,
        use_mean: bool = True,
        use_vmap: Literal["auto", True, False] = "auto",
    ) -> xr.DataArray:
        """Compute the local sample representation of the cells in the ``adata`` object.

        For each cell, it returns a matrix of size ``(n_sample, n_features)``.

        Parameters
        ----------
        adata
            AnnData object to use for computing the local sample representation.
        batch_size
            Batch size to use for computing the local sample representation.
        use_mean
            Whether to use the mean of the latent representation as the local sample
            representation.
        use_vmap
            Whether to use vmap for computing the local sample representation.
            Disabling vmap can be useful if running out of memory on a GPU.
        """
        reductions = [
            MRVIReduction(
                name="sample_representations",
                input="mean_representations" if use_mean else "sampled_representations",
                fn=lambda x: x,
                group_by=None,
            )
        ]
        return self.compute_local_statistics(
            reductions,
            adata=adata,
            indices=indices,
            batch_size=batch_size,
            use_vmap=use_vmap,
        ).sample_representations

    def get_local_sample_distances(
        self,
        adata: AnnData | None = None,
        batch_size: int = 256,
        use_mean: bool = True,
        normalize_distances: bool = False,
        use_vmap: Literal["auto", True, False] = "auto",
        groupby: list[str] | str | None = None,
        keep_cell: bool = True,
        norm: str = "l2",
        mc_samples: int = 10,
    ) -> xr.Dataset:
        """Compute local sample distances.

        Computes cell-specific distances between samples, of size ``(n_sample, n_sample)``,
        stored as a Dataset, with variable name ``"cell"``, of size
        ``(n_cell, n_sample, n_sample)``. If in addition, ``groupby`` is provided, distances are
        also aggregated by group. In this case, the group-specific distances via group name key.

        Parameters
        ----------
        adata
            AnnData object to use for computing the local sample representation.
        batch_size
            Batch size to use for computing the local sample representation.
        use_mean
            Whether to use the mean of the latent representation as the local sample
            representation.
        normalize_distances
            Whether to normalize the local sample distances. Normalizes by the standard deviation
            of the original intra-sample distances. Only works with ``use_mean=False``.
        use_vmap
            Whether to use vmap for computing the local sample representation. Disabling vmap can
            be useful if running out of memory on a GPU.
        groupby
            List of categorical keys or single key of the anndata that is used to group the cells.
        keep_cell
            Whether to keep the original cell sample-sample distance matrices.
        norm
            Norm to use for computing the local sample distances.
        mc_samples
            Number of Monte Carlo samples to use for computing the local sample distances. Only
            relevant if ``use_mean=False``.
        """
        use_vmap = "auto" if use_vmap == "auto" else use_vmap

        input = "mean_distances" if use_mean else "sampled_distances"
        if normalize_distances:
            if use_mean:
                warnings.warn(
                    "Normalizing distances uses sampled distances. Ignoring ``use_mean``.",
                    UserWarning,
                    stacklevel=2,
                )
            input = "normalized_distances"
        if groupby and not isinstance(groupby, list):
            groupby = [groupby]

        reductions = []
        if not keep_cell and not groupby:
            raise ValueError("Undefined computation because not keep_cell and no groupby.")
        if keep_cell:
            reductions.append(
                MRVIReduction(
                    name="cell",
                    input=input,
                    fn=lambda x: x,
                )
            )
        if groupby:
            for groupby_key in groupby:
                reductions.append(
                    MRVIReduction(
                        name=groupby_key,
                        input=input,
                        group_by=groupby_key,
                    )
                )
        return self.compute_local_statistics(
            reductions,
            adata=adata,
            batch_size=batch_size,
            use_vmap=use_vmap,
            norm=norm,
            mc_samples=mc_samples,
        )

    def get_aggregated_posterior(
        self,
        adata: AnnData | None = None,
        sample: str | int | None = None,
        indices: npt.ArrayLike | None = None,
        batch_size: int = 256,
    ) -> Distribution:
        """Computes the aggregated posterior over the ``u`` latent representations.

        For the specified samples, it computes the aggregated posterior over the ``u`` latent
        representations. Returns a NumPyro MixtureSameFamily distribution.

        Parameters
        ----------
        adata
            AnnData object to use. Defaults to the AnnData object used to initialize the model.
        sample
            Name or index of the sample to filter on. If ``None``, uses all cells.
        indices
            Indices of cells to use.
        batch_size
            Batch size to use for computing the latent representation.

        Returns
        -------
        A mixture distribution of the aggregated posterior.
        """
        from torch.distributions import (
            Categorical,
            MixtureSameFamily,
            Normal,
        )

        self._check_if_trained(warn=False)
        adata = self._validate_anndata(adata)
        if indices is None:
            indices = np.arange(adata.n_obs)
        if sample is not None:
            indices = np.intersect1d(
                np.array(indices), np.where(adata.obs[self.sample_key] == sample)[0]
            )
        dataloader = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size, iter_ndarray=True
        )

        qu_locs = []
        qu_scales = []
        for tensor in dataloader:
            inference_inputs = self.module._get_inference_input(tensor)
            inference_inputs.update({"use_mean": True, "sample_index": tensor["ind_x"]})
            use_mean = inference_inputs.get("use_mean", False)

            inference_inputs = self.module._get_inference_input(tensor)
            inference_inputs["use_mean"] = use_mean
            # inference_inputs["sample_index"] = tensor["ind_x"]
            outputs = self.module.inference(**inference_inputs)

            qu_locs.append(outputs["qu"].loc)
            qu_scales.append(outputs["qu"].scale)

        qu_loc = torch.cat(qu_locs, axis=0)  # n_cells x n_latent_u
        qu_scale = torch.cat(qu_scales, axis=0)  # n_cells x n_latent_u
        return MixtureSameFamily(
            Categorical(probs=torch.ones(qu_loc.shape[0]) / qu_loc.shape[0]),
            Independent(Normal(qu_loc, qu_scale), 1),
        )

    def differential_abundance(
        self,
        adata: AnnData | None = None,
        sample_cov_keys: list[str] | None = None,
        sample_subset: list[str] | None = None,
        compute_log_enrichment: bool = False,
        omit_original_sample: bool = True,
        batch_size: int = 128,
    ) -> xr.Dataset:
        """Compute the differential abundance between samples.

        Computes the logarithm of the ratio of the probabilities of each sample conditioned on the
        estimated aggregate posterior distribution of each cell.

        Parameters
        ----------
        adata
            The data object to compute the differential abundance for.
            If not given, the data object stored in the model is used.
        sample_cov_keys
            Keys for covariates (batch, etc.) that should also be taken into account
            when computing the differential abundance. At the moment, only discrete covariates are
            supported.
        sample_subset
            Only compute differential abundance for these sample labels.
        compute_log_enrichment
            Whether to compute the log enrichment scores for each covariate value.
        omit_original_sample
            If true, each cell's sample-of-origin is discarded to compute aggregate posteriors.
            Only relevant if sample_cov_keys is not None.
        batch_size
            Minibatch size for computing the differential abundance.

        Returns
        -------
        A dataset with data variables:

        * ``"log_probs"``: Array of shape ``(n_cells, n_samples)`` containing the log probabilities
            for each cell across samples.
        * ``"{cov_key}_log_probs"``: For each key in ``sample_cov_keys``, an array of shape
            ``(n_cells, _cov_values)`` containing the log probabilities for each cell across
            covariate values.
        """
        from pandas import DataFrame
        from scipy.special import logsumexp

        adata = self._validate_anndata(adata)

        if sample_cov_keys is not None:
            for key in sample_cov_keys:
                n_cov_values = len(adata.obs[key].unique())
                n_samples = len(adata.obs[self.sample_key].unique())
                if n_cov_values > n_samples / 2:
                    warnings.warn(
                        f"The covariate '{key}' does not seem to refer to a discrete key. "
                        f"It has {n_cov_values} unique values, which exceeds one half of the "
                        f"total samples ({n_samples}).",
                        UserWarning,
                        stacklevel=2,
                    )

        us = self.get_latent_representation(
            adata, use_mean=True, give_z=False, batch_size=batch_size
        )

        log_probs = []
        unique_samples = adata.obs[self.sample_key].unique()
        for sample_name in tqdm(unique_samples):
            ap = self.get_aggregated_posterior(
                adata=adata, sample=sample_name, batch_size=batch_size
            )
            n_splits = max(adata.n_obs // batch_size, 1)
            log_probs_ = []
            for u_rep in np.array_split(us, n_splits):
                if not isinstance(u_rep, torch.Tensor):
                    u_rep = torch.tensor(u_rep, dtype=torch.float32, device=self.device)
                log_probs_.append(ap.log_prob(u_rep).detach().cpu().numpy()[..., np.newaxis])

            log_probs.append(np.concatenate(log_probs_, axis=0))  # (n_cells, 1)

        log_probs = np.concatenate(log_probs, 1)

        coords = {
            "cell_name": adata.obs_names.to_numpy(),
            "sample": unique_samples,
        }
        data_vars = {
            "log_probs": (["cell_name", "sample"], log_probs),
        }
        log_probs_arr = xr.Dataset(data_vars, coords=coords)

        if sample_cov_keys is None or len(sample_cov_keys) == 0:
            return log_probs_arr

        def aggregate_log_probs(log_probs, samples, omit_original_sample=False):
            sample_log_probs = log_probs.loc[
                {"sample": samples}
            ].values  # (n_cells, n_samples_in_group)
            if omit_original_sample:
                sample_one_hot = np.zeros((adata.n_obs, len(samples)))
                for i, sample in enumerate(samples):
                    sample_one_hot[adata.obs[self.sample_key] == sample, i] = 1
                log_probs_no_original = np.where(
                    sample_one_hot, -np.inf, sample_log_probs
                )  # virtually discards samples-of-origin from aggregate posteriors
                return logsumexp(log_probs_no_original, axis=1) - np.log(
                    (1 - sample_one_hot).sum(axis=1)
                )
            else:
                return logsumexp(sample_log_probs, axis=1) - np.log(sample_log_probs.shape[1])

        sample_cov_log_probs_map = {}
        sample_cov_log_enrichs_map = {}
        for sample_cov_key in sample_cov_keys:
            sample_cov_unique_values = self.sample_info[sample_cov_key].unique()
            per_val_log_probs = {}
            per_val_log_enrichs = {}
            for sample_cov_value in sample_cov_unique_values:
                cov_samples = (
                    self.sample_info[self.sample_info[sample_cov_key] == sample_cov_value]
                )[self.sample_key].to_numpy()
                if sample_subset is not None:
                    cov_samples = np.intersect1d(cov_samples, np.array(sample_subset))
                if len(cov_samples) == 0:
                    continue

                val_log_probs = aggregate_log_probs(
                    log_probs_arr.log_probs,
                    cov_samples,
                    omit_original_sample=omit_original_sample,
                )
                per_val_log_probs[sample_cov_value] = val_log_probs

                if compute_log_enrichment:
                    rest_samples = np.setdiff1d(unique_samples, cov_samples)
                    if len(rest_samples) == 0:
                        warnings.warn(
                            f"All samples have {sample_cov_key}={sample_cov_value}. Skipping log "
                            "enrichment computation.",
                            UserWarning,
                            stacklevel=2,
                        )
                        continue
                    rest_val_log_probs = aggregate_log_probs(
                        log_probs_arr.log_probs,
                        rest_samples,
                        omit_original_sample=omit_original_sample,
                    )
                    enrichment_scores = val_log_probs - rest_val_log_probs
                    per_val_log_enrichs[sample_cov_value] = enrichment_scores
            sample_cov_log_probs_map[sample_cov_key] = DataFrame.from_dict(per_val_log_probs)
            if compute_log_enrichment and len(per_val_log_enrichs) > 0:
                sample_cov_log_enrichs_map[sample_cov_key] = DataFrame.from_dict(
                    per_val_log_enrichs
                )

        coords = {
            "cell_name": adata.obs_names.to_numpy(),
            "sample": unique_samples,
            **{
                sample_cov_key: sample_cov_log_probs.columns
                for sample_cov_key, sample_cov_log_probs in sample_cov_log_probs_map.items()
            },
        }
        data_vars = {
            "log_probs": (["cell_name", "sample"], log_probs),
            **{
                f"{sample_cov_key}_log_probs": (
                    ["cell_name", sample_cov_key],
                    sample_cov_log_probs.values,
                )
                for sample_cov_key, sample_cov_log_probs in sample_cov_log_probs_map.items()
            },
        }
        if compute_log_enrichment:
            data_vars.update(
                {
                    f"{sample_key}_log_enrichs": (
                        ["cell_name", sample_key],
                        sample_log_enrichs.values,
                    )
                    for sample_key, sample_log_enrichs in sample_cov_log_enrichs_map.items()
                }
            )
        return xr.Dataset(data_vars, coords=coords)

    def get_outlier_cell_sample_pairs(
        self,
        adata: AnnData | None = None,
        subsample_size: int = 5_000,
        quantile_threshold: float = 0.05,
        admissibility_threshold: float = 0.0,
        batch_size: int = 256,
    ) -> xr.Dataset:
        """Compute admissibility scores for cell-sample pairs.

        This function computes the posterior distribution for u for each cell. Then, for every
        cell, it computes the log-probability of the cell under the posterior of each cell
        each sample and takes the maximum value for a given sample as a measure of admissibility
        for that sample. Additionally, it computes a threshold that determines if
        a cell-sample pair is admissible based on the within-sample admissibility scores.

        Parameters
        ----------
        adata
            AnnData object containing the cells for which to compute the outlier cell-sample pairs.
        subsample_size
            Number of cells to use from each sample to approximate the posterior. If None, uses all
            of the available cells.
        quantile_threshold
            Quantile of the within-sample log probabilities to use as a baseline for admissibility.
        admissibility_threshold
            Threshold for admissibility. Cell-sample pairs with admissibility below this threshold
            are considered outliers.
        batch_size
            Size of the batch to use for computing outlier cell-sample pairs.
        """
        adata = self._validate_anndata(adata)
        us = self.get_latent_representation(adata, use_mean=True, give_z=False)
        adata.obsm["U"] = us

        log_probs = []
        threshs = []
        unique_samples = adata.obs[self.sample_key].unique()
        for sample_name in tqdm(unique_samples):
            sample_idxs = np.where(adata.obs[self.sample_key] == sample_name)[0]
            if subsample_size is not None and sample_idxs.shape[0] > subsample_size:
                sample_idxs = np.random.choice(sample_idxs, size=subsample_size, replace=False)
            adata_s = adata[sample_idxs]

            ap = self.get_aggregated_posterior(adata=adata, indices=sample_idxs)

            in_max_comp_log_probs = ap.component_distribution.log_prob(
                torch.unsqueeze(
                    torch.tensor(adata_s.obsm["U"], dtype=torch.float32, device=self.device),
                    dim=-2,
                )
            )
            log_probs_s = rowwise_max_excluding_diagonal(in_max_comp_log_probs)

            log_probs_ = []
            n_splits = adata.n_obs // batch_size
            for u_rep in np.array_split(adata.obsm["U"], n_splits):
                u_rep_tensor = torch.tensor(
                    u_rep,
                    dtype=ap.component_distribution.base_dist.loc.dtype,
                    device=ap.component_distribution.base_dist.loc.device,
                )

                expanded_u_rep = u_rep_tensor.unsqueeze(
                    -2
                )  # Equivalent to np.expand_dims(..., -2)

                log_probs_.append(
                    ap.component_distribution.log_prob(
                        expanded_u_rep  # (n_cells_batch, 1, n_latent_dim)
                    )
                    .max(  # (n_cells_batch, n_cells_ap)
                        axis=1, keepdims=True
                    )
                    .values.detach()
                    .cpu()
                    .numpy()  # (n_cells_batch, 1)
                )

            log_probs_ = np.concatenate(log_probs_, axis=0)  # (n_cells, 1)

            threshs.append(np.array(log_probs_s.detach()))
            log_probs.append(np.array(log_probs_))

        threshs_all = np.concatenate(threshs)
        global_thresh = np.quantile(threshs_all, q=quantile_threshold)
        threshs = np.array(len(log_probs) * [global_thresh])

        log_probs = np.concatenate(log_probs, 1)
        log_ratios = log_probs - threshs

        coords = {
            "cell_name": adata.obs_names.to_numpy(),
            "sample": unique_samples,
        }
        data_vars = {
            "log_probs": (["cell_name", "sample"], log_probs),
            "log_ratios": (
                ["cell_name", "sample"],
                log_ratios,
            ),
            "is_admissible": (
                ["cell_name", "sample"],
                log_ratios > admissibility_threshold,
            ),
        }
        return xr.Dataset(data_vars, coords=coords)

    def differential_expression(
        self,
        adata: AnnData | None = None,
        sample_cov_keys: list[str] | None = None,
        sample_subset: list[str] | None = None,
        batch_size: int = 128,
        use_vmap: Literal["auto", True, False] = "auto",
        normalize_design_matrix: bool = True,
        add_batch_specific_offsets: bool = False,
        mc_samples: int = 100,
        store_lfc: bool = False,
        store_lfc_metadata_subset: list[str] | None = None,
        store_baseline: bool = False,
        eps_lfc: float = 1e-4,
        filter_inadmissible_samples: bool = False,
        lambd: float = 0.0,  # TODO: should it be 0?
        delta: float | None = 0.3,
        **filter_samples_kwargs,
    ) -> xr.Dataset:
        """Compute cell-specific multivariate differential expression.

        For every cell, we first compute all counterfactual cell-state shifts, defined as
        ``e_d = z_d - u``, where ``z_d`` is the latent representation of the cell for sample ``d``
        and ``u`` is the sample-unaware latent representation. Then, we fit a linear model in each
        cell of the form: ``e_d = X_d * beta_d + iid gaussian noise``.

        Parameters
        ----------
        sample_cov_keys
            List of sample covariates to consider for the multivariate analysis.
            These keys should be present in ``adata.obs``.
        adata
            AnnData object to use for computing the local sample representation.
            If ``None``, the analysis is performed on all cells in the dataset.
        sample_subset
            Optional list of samples to consider for the multivariate analysis.
            If ``None``, all samples are considered.
        batch_size
            Batch size to use for computing the local sample representation.
        use_vmap
            Whether to use vmap for computing the local sample representation.
        normalize_design_matrix
            Whether to normalize the design matrix.
        add_batch_specific_offsets
            Whether to offset the design matrix by adding batch-specific offsets to the design
            matrix. Setting this option to True is recommended when considering multi-site
            datasets.
        mc_samples
            How many MC samples should be taken for computing betas.
        store_lfc
            Whether to store the log-fold changes in the module.
            Storing log-fold changes is memory-intensive and may require to specify
            a smaller set of cells to analyze, e.g., by specifying ``adata``.
        store_lfc_metadata_subset
            Specifies a subset of metadata for which log-fold changes are computed.
            These keys must be a subset of ``sample_cov_keys``.
            Only applies when ``store_lfc=True``.
        store_baseline
            Whether to store the expression in the module if logfoldchanges are computed.
        eps_lfc
            Epsilon to add to the log-fold changes to avoid detecting genes with low expression.
        filter_inadmissible_samples
            Whether to filter out-of-distribution samples prior to performing the analysis.
        lambd
            Regularization parameter for the linear model.
        delta
            LFC threshold used to compute posterior DE probabilities.
            If None does not compute them to save memory consumption.
        filter_samples_kwargs
            Keyword arguments to pass to :meth:`~scvi.external.MRVI.get_outlier_cell_sample_pairs`.

        Returns
        -------
        A dataset containing the results of the differential expression analysis:

        * ``"beta"``: Coefficients for each covariate across cells and latent dimensions.
        * ``"effect_size"``: Effect sizes for each covariate across cells.
        * ``"pvalue"``: P-values for each covariate across cells.
        * ``"padj"``: Adjusted P-values for each covariate across cells using the
            Benjamini-Hochberg procedure.
        * ``"lfc"``: Log fold changes for each covariate across cells and genes, if ``store_lfc``
            is ``True``.
        * ``"lfc_std"``: Standard deviation of log fold changes, if ``store_lfc`` is ``True`` and
            ``delta`` is not ``None``.
        * ``"pde"``: Posterior DE probabilities, if ``store_lfc`` is ``True`` and ``delta`` is not
            ``None``.
        * ``"baseline_expression"``: Baseline expression levels for each covariate across cells and
            genes, if ``store_baseline`` is ``True``.
        * ``"n_samples"``: Number of admissible samples for each cell, if
            ``filter_inadmissible_samples`` is ``True``.
        """
        from functools import partial

        import scipy.linalg
        import torch.distributions as dist
        from scipy.stats import false_discovery_control

        use_vmap = use_vmap if use_vmap != "auto" else self.summary_stats.n_sample < 500

        if sample_cov_keys is None:
            # Hack: kept as kwarg to maintain order of arguments.
            raise ValueError("Must assign `sample_cov_keys`")
        adata = self.adata if adata is None else adata
        self._check_if_trained(warn=False)
        # Hack to ensure new AnnDatas have indices and indices have correct dimensions.
        if adata is not None:
            adata.obs["_indices"] = np.arange(adata.n_obs).astype(int)

        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(
            adata=adata,
            indices=None,
            batch_size=batch_size,
            iter_ndarray=True,
        )
        n_sample = self.summary_stats.n_sample
        # vars_in = {"params": self.module.params, **self.module.state}

        sample_mask = (
            np.isin(self.sample_order, sample_subset)
            if sample_subset is not None
            else np.ones(n_sample, dtype=bool)
        )
        sample_mask = np.array(sample_mask)
        sample_order = self.sample_order[sample_mask]
        n_samples_kept = sample_mask.sum()

        if filter_inadmissible_samples:
            admissible_samples = self.get_outlier_cell_sample_pairs(
                adata=adata, **filter_samples_kwargs
            )["is_admissible"].loc[{"sample": sample_order}]
            assert (admissible_samples.sample == sample_order).all()
            admissible_samples = admissible_samples.values
        else:
            admissible_samples = np.ones((adata.n_obs, n_samples_kept), dtype=bool)
        n_admissible_samples = admissible_samples.sum(1)

        (
            Xmat,
            Xmat_names,
            covariates_require_lfc,
            offset_indices,
        ) = self._construct_design_matrix(
            sample_cov_keys=sample_cov_keys,
            sample_mask=sample_mask,
            normalize_design_matrix=normalize_design_matrix,
            add_batch_specific_offsets=add_batch_specific_offsets,
            store_lfc=store_lfc,
            store_lfc_metadata_subset=store_lfc_metadata_subset,
        )
        add_batch_specific_offsets = offset_indices is not None
        n_covariates = Xmat.shape[1]

        def sqrtm_batch(xtmx: torch.Tensor) -> torch.Tensor:
            # xtmx shape: (n_cells, n_covariates, n_covariates)
            sqrtm_list = [
                torch.from_numpy(scipy.linalg.sqrtm(m.cpu().numpy()).real).to(xtmx.device)
                for m in xtmx
            ]
            return torch.stack(sqrtm_list, dim=0)

        def process_design_matrix(
            admissible_samples_dmat: torch.Tensor,
            Xmat: torch.Tensor,
        ) -> tuple[torch.Tensor, torch.Tensor]:
            xtmx = torch.einsum("ak,nkl,lm->nam", Xmat.T, admissible_samples_dmat, Xmat)
            xtmx += lambd * torch.eye(n_covariates)

            prefactor = sqrtm_batch(xtmx)

            inv_ = torch.vmap(torch.linalg.pinv)(xtmx)
            Amat = torch.einsum("nab,bc,ncd->nad", inv_, Xmat.T, admissible_samples_dmat)
            return Amat, prefactor

        def mapped_inference_fn(
            x: torch.Tensor,
            sample_index: torch.Tensor,
            cf_sample: torch.Tensor,
            Amat: torch.Tensor,
            prefactor: torch.Tensor,
            n_samples_per_cell: int,
            admissible_samples_mat: torch.Tensor,
            use_mean: bool,
            mc_samples: int,
        ):
            def inference_fn(cf):
                return self.module.inference(
                    x, sample_index, cf_sample=cf, mc_samples=mc_samples, use_mean=use_mean
                )["z"]

            inference_partial = partial(
                self.module.inference,
                x=x,
                sample_index=sample_index,
                mc_samples=mc_samples,
                use_mean=use_mean,
            )

            if use_vmap:
                eps_ = torch.vmap(inference_fn, in_dims=0, out_dims=-2, randomness="different")(
                    cf_sample
                )
            else:

                def per_sample_inference_fn(cf_sample):
                    return inference_partial(cf_sample=cf_sample)["z"]

                eps_ = torch.stack(
                    [per_sample_inference_fn(cf_sample=cf_sample_i) for cf_sample_i in cf_sample]
                ).permute(1, 2, 0, 3)

            eps_std = eps_.std(axis=2, keepdims=True)
            eps_mean = eps_.mean(axis=2, keepdims=True)

            eps = (eps_ - eps_mean) / (1e-6 + eps_std)  # over samples
            # MLE for betas
            betas = torch.einsum("nks,ansd->ankd", Amat, eps)

            # Statistical tests

            betas_norm = torch.einsum("ankd,nkl->anld", betas, prefactor)
            ts = (betas_norm**2).mean(axis=0).sum(axis=-1)

            chi2_dist = dist.Chi2(n_samples_per_cell[:, None])
            pvals = 1 - chi2_dist.cdf(ts)

            betas = betas * eps_std

            lfc_mean = None
            lfc_std = None
            pde = None
            if store_lfc:
                betas_ = betas.permute((0, 2, 1, 3))
                eps_mean_ = eps_mean.permute((0, 2, 1, 3))
                betas_covariates = betas_[:, covariates_require_lfc, :, :]

                def h_inference_fn(extra_eps, batch_index_cf, batch_offset_eps):
                    extra_eps = extra_eps + batch_offset_eps

                    return self.module.compute_h_from_x_eps(
                        x=x,
                        sample_index=sample_index,
                        extra_eps=extra_eps,
                        batch_index=batch_index_cf,
                        cf_sample=None,
                        mc_samples=None,
                        # mc_samples also taken for eps. vmap over mc_samples
                    )

                batch_index_ = torch.arange(self.summary_stats.n_batch)[:, None]
                batch_index_ = batch_index_.repeat(1, n_cells)[..., None]  # (n_batch, n_cells, 1)
                betas_null = torch.zeros_like(betas_covariates)

                if add_batch_specific_offsets:
                    batch_weights = torch.einsum(
                        "nd,db->nb", admissible_samples_mat, Xmat[:, offset_indices]
                    ).mean(0)
                    betas_offset_ = betas_[:, offset_indices, :, :] + eps_mean_
                else:
                    batch_weights = (1.0 / self.summary_stats.n_batch) * torch.ones(
                        self.summary_stats.n_batch
                    )
                    mc_samples, _, n_cells_, n_latent = betas_covariates.shape
                    betas_offset_ = (
                        torch.zeros((mc_samples, self.summary_stats.n_batch, n_cells_, n_latent))
                        + eps_mean_
                    )
                # batch_offset shape (mc_samples, n_batch, n_cells, n_latent)

                f_ = torch.vmap(
                    h_inference_fn, in_dims=(0, None, 0), out_dims=0, randomness="different"
                )  # fn over MC samples
                f_ = torch.vmap(
                    f_, in_dims=(1, None, None), out_dims=1, randomness="different"
                )  # fn over covariates
                f_ = torch.vmap(
                    f_, in_dims=(None, 0, 1), out_dims=0, randomness="different"
                )  # fn over batches
                h_fn = f_

                x_1 = h_fn(betas_covariates, batch_index_, betas_offset_)
                x_0 = h_fn(betas_null, batch_index_, betas_offset_)

                lfcs = torch.log2(x_1 + eps_lfc) - torch.log2(x_0 + eps_lfc)
                # Compute weighted average manually: sum(weights * values) / sum(weights)
                lfcs_mean_over_mc = lfcs.mean(1)  # (n_batch, n_covariates, n_cells, n_genes)
                lfc_mean = (batch_weights[:, None, None, None] * lfcs_mean_over_mc).sum(
                    0
                ) / batch_weights.sum()
                if delta is not None:
                    lfcs_var_over_mc = lfcs.var(1)  # (n_batch, n_covariates, n_cells, n_genes)
                    lfc_std = torch.sqrt(
                        (batch_weights[:, None, None, None] * lfcs_var_over_mc).sum(0)
                        / batch_weights.sum()
                    )
                    pde = (torch.abs(lfcs) >= delta).float().mean(1).mean(0)

            if store_baseline:
                baseline_expression = x_1.mean(1)
            else:
                baseline_expression = None
            return {
                "beta": betas.mean(0),
                "effect_size": ts,
                "pvalue": pvals,
                "lfc_mean": lfc_mean,
                "lfc_std": lfc_std,
                "pde": pde,
                "baseline_expression": baseline_expression,
            }

        beta = []
        effect_size = []
        pvalue = []
        lfc = []
        lfc_std = []
        pde = []
        baseline_expression = []
        for tensors in tqdm(scdl):
            indices = tensors[REGISTRY_KEYS.INDICES_KEY].astype(int).flatten()
            n_cells = tensors[REGISTRY_KEYS.X_KEY].shape[0]
            cf_sample = np.broadcast_to(
                (np.where(sample_mask)[0])[:, None, None], (n_samples_kept, n_cells, 1)
            )
            inf_inputs = self.module._get_inference_input(
                tensors,
            )
            admissible_samples_mat = torch.Tensor(
                admissible_samples[indices]
            )  # (n_cells, n_samples)
            n_samples_per_cell = admissible_samples_mat.sum(axis=1)
            admissible_samples_dmat = torch.vmap(torch.diag)(
                admissible_samples_mat
            ).float()  # (n_cells, n_samples, n_samples)
            # element nij is 1 if sample i is admissible and i=j for cell n
            Amat, prefactor = process_design_matrix(admissible_samples_dmat, Xmat)
            Amat = Amat.to(self.device)
            prefactor = prefactor.to(self.device)

            try:
                res = mapped_inference_fn(
                    x=torch.Tensor(inf_inputs["x"]),
                    sample_index=torch.Tensor(inf_inputs["sample_index"]),
                    cf_sample=torch.Tensor(cf_sample),
                    Amat=Amat,
                    prefactor=prefactor,
                    n_samples_per_cell=n_samples_per_cell,
                    admissible_samples_mat=admissible_samples_mat,
                    use_mean=False,
                    mc_samples=mc_samples,
                )
            except RuntimeError as e:
                if use_vmap:
                    raise RuntimeError(
                        "Out of memory. Try setting use_vmap=False."
                    ) from e  # TODO: update error msg
                else:
                    raise e

            beta.append(np.array(res["beta"].detach().cpu()))
            effect_size.append(np.array(res["effect_size"].detach().cpu()))
            pvalue.append(np.array(res["pvalue"].detach().cpu()))
            if store_lfc:
                lfc.append(
                    np.array(res["lfc_mean"].detach().cpu()).transpose(1, 0, 2)
                )  # Transpose to (n_cells, n_covariates, n_genes)
                if delta is not None:
                    lfc_std.append(
                        np.array(res["lfc_std"].detach().cpu()).transpose(1, 0, 2)
                    )  # Transpose to (n_cells, n_covariates, n_genes)
                    pde.append(
                        np.array(res["pde"].detach().cpu()).transpose(1, 0, 2)
                    )  # Transpose to (n_cells, n_covariates, n_genes)
            if store_baseline:
                baseline_expression.append(np.array(res["baseline_expression"].detach().cpu()))
        beta = np.concatenate(beta, axis=0)  # concatenate along cell dimension
        effect_size = np.concatenate(effect_size, axis=0)  # concatenate along cell dimension
        pvalue = np.concatenate(pvalue, axis=0)  # concatenate along cell dimension
        pvalue_shape = pvalue.shape
        padj = false_discovery_control(pvalue.flatten(), method="bh").reshape(pvalue_shape)

        coords = {
            "cell_name": (("cell_name"), adata.obs_names),
            "covariate": (("covariate"), Xmat_names),
            "latent_dim": (("latent_dim"), np.arange(beta.shape[2])),
            "gene": (("gene"), adata.var_names),
        }
        data_vars = {
            "beta": (
                ["cell_name", "covariate", "latent_dim"],
                beta,
            ),
            "effect_size": (
                ["cell_name", "covariate"],
                effect_size,
            ),
            "pvalue": (
                ["cell_name", "covariate"],
                pvalue,
            ),
            "padj": (
                ["cell_name", "covariate"],
                padj,
            ),
        }
        if filter_inadmissible_samples:
            data_vars["n_samples"] = (
                ["cell_name"],
                n_admissible_samples,
            )
        if store_lfc:
            if store_lfc_metadata_subset is None and not add_batch_specific_offsets:
                coords_lfc = ["cell_name", "covariate", "gene"]
            else:
                coords_lfc = ["cell_name", "covariate_sub", "gene"]
                coords["covariate_sub"] = (
                    ("covariate_sub"),
                    Xmat_names[covariates_require_lfc],
                )
            lfc = np.concatenate(lfc, axis=0)  # concatenate along cell dimension
            data_vars["lfc"] = (coords_lfc, lfc)
            if delta is not None:
                lfc_std = np.concatenate(lfc_std, axis=0)  # concatenate along cell dimension
                pde = np.concatenate(pde, axis=0)  # concatenate along cell dimension
                data_vars["lfc_std"] = (coords_lfc, lfc_std)
                data_vars["pde"] = (coords_lfc, pde)

        if store_baseline:
            baseline_expression = np.concatenate(baseline_expression, axis=1)
            data_vars["baseline_expression"] = (
                ["covariate", "cell_name", "gene"],
                baseline_expression,
            )
        return xr.Dataset(data_vars, coords=coords)

    def _construct_design_matrix(
        self,
        sample_cov_keys: list[str],
        sample_mask: npt.ArrayLike,
        normalize_design_matrix: bool,
        add_batch_specific_offsets: bool,
        store_lfc: bool,
        store_lfc_metadata_subset: list[str] | None = None,
    ) -> tuple[torch.Tensor, npt.NDArray, torch.Tensor, torch.Tensor | None]:
        """Construct a design matrix of samples and covariates.

        Starting from a list of sample covariate keys, construct a design matrix of samples and
        covariates. Categorical covariates are one-hot encoded.

        Parameters
        ----------
        sample_cov_keys
            List of sample metadata to use as covariates.
        sample_mask
            Mask of admissible samples. Must have the same length as the number of samples in the
            dataset.
        normalize_design_matrix
            Whether the design matrix should be 0-1 normalized. This is useful to ensure that the
            beta coefficients are comparable across covariates.
        add_batch_specific_offsets
            Whether the design matrix should be offset. If True, the matrix includes batch-specific
            offsets. This ensures that we can learn perturbation effects that do not depend on
            batch effects.

        Returns
        -------
        A tuple consisting of:

        1. The design matrix
        2. Names for each column in the design matrix
        3. A mask precising which coefficients from the design matrix require to compute LFCs.
        4. A mask precising which coefficients from the design matrix correspond to offsets.
        """
        from pandas import Series, get_dummies

        Xmat = []
        Xmat_names = []
        Xmat_dim_to_key = []
        sample_info = self.sample_info.iloc[sample_mask]
        for sample_cov_key in tqdm(sample_cov_keys):
            cov = sample_info[sample_cov_key]
            if (cov.dtype == str) or (cov.dtype == "category"):
                cov = cov.cat.remove_unused_categories()
                cov = get_dummies(cov, drop_first=True)
                cov_names = np.array([f"{sample_cov_key}_{col}" for col in cov.columns])
                cov = cov.values
            else:
                cov_names = np.array([sample_cov_key])
                cov = cov.values[:, None]
            n_covs = cov.shape[1]
            Xmat.append(cov)
            Xmat_names.append(cov_names)
            Xmat_dim_to_key.append([sample_cov_key] * n_covs)
        Xmat_names = np.concatenate(Xmat_names)
        Xmat = np.concatenate(Xmat, axis=1).astype(np.float32)
        Xmat_dim_to_key = np.concatenate(Xmat_dim_to_key)

        if normalize_design_matrix:
            Xmat = (Xmat - Xmat.min(axis=0)) / (1e-6 + Xmat.max(axis=0) - Xmat.min(axis=0))
        if add_batch_specific_offsets:
            cov = sample_info["_scvi_batch"]
            if cov.nunique() == self.summary_stats.n_batch:
                cov = np.eye(self.summary_stats.n_batch)[sample_info["_scvi_batch"].values]
                cov_names = ["offset_batch_" + str(i) for i in range(self.summary_stats.n_batch)]
                Xmat = np.concatenate([cov, Xmat], axis=1)
                Xmat_names = np.concatenate([np.array(cov_names), Xmat_names])
                Xmat_dim_to_key = np.concatenate([np.array(cov_names), Xmat_dim_to_key])

                # Retrieve indices of offset covariates in the right order
                offset_indices = (
                    Series(np.arange(len(Xmat_names)), index=Xmat_names).loc[cov_names].values
                )
                offset_indices = torch.Tensor(offset_indices)
            else:
                warnings.warn(
                    """
                        Number of batches in sample_info does not match number of batches in
                        summary_stats. `add_batch_specific_offsets=True` assumes that samples are
                        not shared across batches. Setting `add_batch_specific_offsets=False`...
                    """,
                    stacklevel=2,
                )
                offset_indices = None
        else:
            offset_indices = None

        Xmat = torch.Tensor(Xmat)
        if store_lfc:
            covariates_require_lfc = (
                np.isin(Xmat_dim_to_key, store_lfc_metadata_subset)
                if store_lfc_metadata_subset is not None
                else np.isin(Xmat_dim_to_key, sample_cov_keys)
            )
        else:
            covariates_require_lfc = np.zeros(len(Xmat_names), dtype=bool)
        covariates_require_lfc = torch.Tensor(covariates_require_lfc).bool()
        return Xmat, Xmat_names, covariates_require_lfc, offset_indices

    def update_sample_info(self, adata):
        """Initialize/update metadata in the case where additional covariates are added.

        Parameters
        ----------
        adata
            AnnData object to update the sample info with. Typically, this corresponds to the
            working dataset, where additional sample-specific covariates have been added.

        Examples
        --------
        >>> import scanpy as sc
        >>> from scvi.external import MRVI
        >>> MRVI.setup_anndata(adata, sample_key="sample_id")
        >>> model = MRVI(adata)
        >>> model.train()
        >>> # Update sample info with new covariates
        >>> sample_mapper = {"sample_1": "healthy", "sample_2": "disease"}
        >>> adata.obs["disease_status"] = adata.obs["sample_id"].map(sample_mapper)
        >>> model.update_sample_info(adata)
        """
        adata = self._validate_anndata(adata)
        obs_df = adata.obs.copy()
        obs_df = obs_df.loc[~obs_df._scvi_sample.duplicated("first")]
        self.sample_info = obs_df.set_index("_scvi_sample").sort_index()
