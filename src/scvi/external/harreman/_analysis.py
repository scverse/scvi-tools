from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from scvi.external.harreman._constants import (
    HARREMAN_DENOISED_LAYER,
    HARREMAN_LATENT_OBSM,
    HARREMAN_PARAMS_KEY,
    HARREMAN_UNS_KEY,
    MODEL_DESTVI,
    MODEL_RESOLVI,
    MODEL_SCVIVA,
    STEP_CCC,
    STEP_FILTER,
    STEP_GENE_PAIRS,
    STEP_ICS,
    STEP_SETUP,
    STEP_SIG,
    SUPPORTED_MODELS,
)
from scvi.external.harreman._results import HarremanResults
from scvi.external.harreman.preprocessing.database import (
    extract_interaction_db as _extract_interaction_db,
)
from scvi.external.harreman.tools.cell_communication import (
    apply_gene_filtering as _apply_gene_filtering,
)
from scvi.external.harreman.tools.cell_communication import (
    compute_cell_communication as _compute_cell_communication,
)
from scvi.external.harreman.tools.cell_communication import (
    compute_ct_cell_communication as _compute_ct_cell_communication,
)
from scvi.external.harreman.tools.cell_communication import (
    compute_ct_interacting_cell_scores as _compute_ct_interacting_cell_scores,
)
from scvi.external.harreman.tools.cell_communication import (
    compute_gene_pairs as _compute_gene_pairs,
)
from scvi.external.harreman.tools.cell_communication import (
    compute_interacting_cell_scores as _compute_interacting_cell_scores,
)
from scvi.external.harreman.tools.cell_communication import (
    select_significant_interactions as _select_significant_interactions,
)
from scvi.external.harreman.tools.knn import compute_knn_graph as _compute_knn_graph

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData

    from scvi.model.base import BaseModelClass

logger = logging.getLogger(__name__)

_STEP_PREREQUISITES: dict[str, str | None] = {
    STEP_SETUP: None,
    STEP_FILTER: STEP_SETUP,
    STEP_GENE_PAIRS: STEP_SETUP,
    STEP_CCC: STEP_GENE_PAIRS,
    STEP_ICS: STEP_CCC,
    STEP_SIG: STEP_CCC,
}


class HarremanAnalysis:
    """Downstream spatial metabolic cell-cell communication analysis.

    Parameters
    ----------
    adata
        Spatial AnnData. Must have neighbor coordinates in ``obsm``.
    model
        Optional trained scvi spatial model (RESOLVI, SCVIVA, DestVI).

        - DestVI: calls ``model.get_proportions()`` and attaches each cell-type
          proportion-weighted count matrix as ``adata.layers[cell_type]``.
        - RESOLVI: calls ``model.get_normalized_expression()`` and attaches as
          ``adata.layers[HARREMAN_DENOISED_LAYER]``.
        - SCVIVA: calls ``model.get_latent_representation()`` and attaches as
          ``adata.obsm[HARREMAN_LATENT_OBSM]``.

        Model reference is dropped after extraction.
    layer_key
        Layer to use for counts. ``None`` uses ``adata.X``.

    Examples
    --------
    Without model:

    >>> ha = HarremanAnalysis(adata)
    >>> ha.setup(compute_neighbors_on_key="spatial", species="human")
    >>> ha.filter_genes()
    >>> ha.compute_gene_pairs()
    >>> ha.compute_cell_communication()
    >>> ha.select_significant_interactions()
    >>> results = ha.results

    With DestVI model:

    >>> ha = HarremanAnalysis(adata, model=destvi_model)
    >>> ha.setup(compute_neighbors_on_key="spatial", cell_type_key="cell_type")
    >>> ha.compute_gene_pairs()
    >>> ha.compute_cell_communication(mode="cell_type")
    """

    def __init__(
        self,
        adata: AnnData,
        model: BaseModelClass | None = None,
        layer_key: str | None = None,
    ) -> None:
        from anndata import AnnData as _AnnData

        if not isinstance(adata, _AnnData):
            raise TypeError(f"Expected AnnData, got {type(adata).__name__}")

        self._adata = adata
        self._layer_key = layer_key
        self._completed_steps: set[str] = set()
        self._is_deconvolved = False
        self._ccc_mode: str | None = None

        if model is not None:
            self._apply_model_integration(model, adata)

        self._adata.uns.setdefault(HARREMAN_UNS_KEY, {HARREMAN_PARAMS_KEY: {}})

        self.hs = _HarremanHsAccessor(self)
        self.tl = _HarremanTlAccessor(self)
        self.pl = _HarremanPlAccessor(self)

    # ── Model integration ──────────────────────────────────────────────────────

    def _apply_model_integration(self, model: BaseModelClass, adata: AnnData) -> None:
        model_type = type(model).__name__
        if model_type not in SUPPORTED_MODELS:
            raise ValueError(
                f"Unsupported model type '{model_type}'. Supported: {SUPPORTED_MODELS}"
            )
        if model_type == MODEL_DESTVI:
            self._extract_destvi_outputs(model, adata)
            self._is_deconvolved = True
        elif model_type == MODEL_RESOLVI:
            self._extract_resolvi_outputs(model, adata)
        elif model_type == MODEL_SCVIVA:
            self._extract_scviva_outputs(model, adata)

    @staticmethod
    def _extract_destvi_outputs(model: BaseModelClass, adata: AnnData) -> None:
        """Attach cell-type proportion-scaled count layers from DestVI."""
        import numpy as np
        from scipy.sparse import issparse

        proportions = model.get_proportions()  # DataFrame (n_obs, n_labels)
        X = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
        for ct in proportions.columns:
            weights = proportions[ct].values[:, None]  # (n_obs, 1)
            adata.layers[ct] = X * weights
        logger.info(
            "HarremanAnalysis: attached %d DestVI cell-type layers.", len(proportions.columns)
        )

    @staticmethod
    def _extract_resolvi_outputs(model: BaseModelClass, adata: AnnData) -> None:
        """Attach denoised counts from RESOLVI as a layer."""
        import numpy as np
        import pandas as pd

        denoised = model.get_normalized_expression()  # DataFrame (n_obs, n_vars)
        if isinstance(denoised, pd.DataFrame):
            adata.layers[HARREMAN_DENOISED_LAYER] = denoised.values
        else:
            adata.layers[HARREMAN_DENOISED_LAYER] = np.asarray(denoised)
        logger.info("HarremanAnalysis: attached RESOLVI denoised layer.")

    @staticmethod
    def _extract_scviva_outputs(model: BaseModelClass, adata: AnnData) -> None:
        """Attach latent representation from SCVIVA for KNN construction."""
        latent = model.get_latent_representation()  # ndarray (n_obs, n_latent)
        adata.obsm[HARREMAN_LATENT_OBSM] = latent
        logger.info("HarremanAnalysis: attached SCVIVA latent representation.")

    # ── Setup ──────────────────────────────────────────────────────────────────

    def setup(
        self,
        compute_neighbors_on_key: str,
        species: Literal["human", "mouse"] = "human",
        database: Literal["transporter", "LR", "both"] = "both",
        n_neighbors: int | None = 30,
        neighborhood_radius: float | None = None,
        cell_type_key: str | None = None,
        sample_key: str | None = None,
        spot_diameter: int = 10,
        extracellular_only: bool = True,
        use_gpu: bool | None = None,
    ) -> None:
        """Register adata with HarremanDB and build the spatial proximity graph.

        Parameters
        ----------
        compute_neighbors_on_key
            Key in ``adata.obsm`` to use for KNN computation (e.g. ``"spatial"``).
        species
            ``"human"`` or ``"mouse"``.
        database
            Which database to load: ``"transporter"``, ``"LR"``, or ``"both"``.
        n_neighbors
            Number of nearest neighbors. Either this or ``neighborhood_radius`` required.
        neighborhood_radius
            Radius for neighborhood graph. Alternative to ``n_neighbors``.
        cell_type_key
            Key in ``adata.obs`` for cell types. Required for ``mode="cell_type"``.
        sample_key
            Key in ``adata.obs`` for sample/batch.
        spot_diameter
            Spot diameter of the spatial technology.
        extracellular_only
            Restrict HarremanDB to extracellular metabolites only.
        use_gpu
            ``None`` (default) tries GPU via rapids_singlecell/cuML and falls back to CPU.
            ``False`` forces CPU. ``True`` raises if GPU is unavailable.
        """
        _extract_interaction_db(
            self._adata,
            species=species,
            database=database,
            extracellular_only=extracellular_only,
        )
        _compute_knn_graph(
            self._adata,
            compute_neighbors_on_key=compute_neighbors_on_key,
            n_neighbors=n_neighbors,
            neighborhood_radius=neighborhood_radius,
            sample_key=sample_key,
            use_gpu=use_gpu,
        )
        # Invalidate any downstream steps so re-running setup with new params doesn't
        # leave stale gene-pair or CCC results silently accepted.
        downstream = {STEP_FILTER, STEP_GENE_PAIRS, STEP_CCC, STEP_ICS, STEP_SIG}
        self._completed_steps -= downstream
        self._ccc_mode = None

        # Store run-specific keys under the harreman namespace to avoid collisions
        # with other scvi tools that write to the same top-level uns keys.
        harreman_ns = self._adata.uns[HARREMAN_UNS_KEY]
        if cell_type_key is not None:
            harreman_ns["cell_type_key"] = cell_type_key
        if sample_key is not None:
            harreman_ns["sample_key"] = sample_key
        harreman_ns["spot_diameter"] = spot_diameter

        harreman_ns[HARREMAN_PARAMS_KEY].update(
            {
                "species": species,
                "database": database,
                "layer_key": self._layer_key,
                "compute_neighbors_on_key": compute_neighbors_on_key,
                "cell_type_key": cell_type_key,
                "sample_key": sample_key,
                "spot_diameter": spot_diameter,
            }
        )
        self._completed_steps.add(STEP_SETUP)
        logger.info(
            "HarremanAnalysis: setup complete (species=%s, database=%s).", species, database
        )

    # ── Prerequisite checking ──────────────────────────────────────────────────

    def _require(self, step: str) -> None:
        prereq = _STEP_PREREQUISITES.get(step)
        if prereq is not None and prereq not in self._completed_steps:
            raise RuntimeError(
                f"{step}() requires {prereq}() to be run first. "
                f"Call ha.{prereq}() before continuing."
            )

    # ── Properties ────────────────────────────────────────────────────────────

    @property
    def adata(self) -> AnnData:
        """The working AnnData object."""
        return self._adata

    @property
    def results(self) -> HarremanResults:
        """Typed view over ``adata.uns['harreman']``."""
        if STEP_SETUP not in self._completed_steps:
            raise RuntimeError("No results available. Run ha.setup() first, then analysis steps.")
        return HarremanResults.from_uns(self._adata.uns[HARREMAN_UNS_KEY])

    @property
    def is_deconvolved(self) -> bool:
        """True if model-extracted cell-type proportion layers are present."""
        return self._is_deconvolved

    @property
    def is_set_up(self) -> bool:
        """True if setup() has been called."""
        return STEP_SETUP in self._completed_steps

    def __repr__(self) -> str:
        steps = sorted(self._completed_steps)
        return (
            f"HarremanAnalysis("
            f"is_set_up={self.is_set_up}, "
            f"is_deconvolved={self.is_deconvolved}, "
            f"completed_steps={steps})"
        )

    # ── Analysis steps ────────────────────────────────────────────────────────

    def filter_genes(
        self,
        model: Literal["danb", "bernoulli", "normal", "none"] = "danb",
        feature_elimination: bool = False,
        threshold: float = 0.2,
        autocorrelation_filt: bool = False,
        expression_filt: bool = False,
        de_filt: bool = False,
        device: str | None = None,
        verbose: bool = False,
    ) -> None:
        """Filter genes before gene pair computation.

        Parameters
        ----------
        model
            Expression model for autocorrelation-based filtering.
        feature_elimination
            Filter genes expressed in fewer than ``threshold`` fraction of cells.
        threshold
            Minimum fraction of cells expressing a gene.
        autocorrelation_filt
            Keep only genes with significant local autocorrelation.
        expression_filt
            Keep only genes expressed in each cell type.
        de_filt
            Keep only genes differentially expressed between cell types.
        device
            PyTorch device string. ``None`` uses auto-selection.
        verbose
            Print progress.
        """
        self._require(STEP_FILTER)
        _apply_gene_filtering(
            adata=self._adata,
            layer_key=self._layer_key,
            cell_type_key=self._adata.uns[HARREMAN_UNS_KEY].get("cell_type_key"),
            model=model,
            feature_elimination=feature_elimination,
            threshold=threshold,
            autocorrelation_filt=autocorrelation_filt,
            expression_filt=expression_filt,
            de_filt=de_filt,
            umi_counts_obs_key=None,
            device=device if device is not None else "auto",
            verbose=verbose,
        )
        self._completed_steps.add(STEP_FILTER)

    def compute_gene_pairs(
        self,
        cell_type_pairs: list | None = None,
        ct_specific: bool = True,
        fix_ct: str | None = None,
        verbose: bool = False,
    ) -> None:
        """Compute local correlations for gene pairs in the metabolite database.

        Requires ``setup()`` to have been called.

        Parameters
        ----------
        cell_type_pairs
            List of ``(ct_a, ct_b)`` tuples. If ``None``, uses all pairs.
        ct_specific
            Restrict gene pairs to cell-type-relevant combinations.
        fix_ct
            Fix one cell type across all pairs.
        verbose
            Print progress messages.
        """
        self._require(STEP_GENE_PAIRS)
        _compute_gene_pairs(
            adata=self._adata,
            layer_key=self._layer_key,
            cell_type_key=self._adata.uns[HARREMAN_UNS_KEY].get("cell_type_key"),
            cell_type_pairs=cell_type_pairs,
            ct_specific=ct_specific,
            fix_ct=fix_ct,
            verbose=verbose,
        )
        self._completed_steps.add(STEP_GENE_PAIRS)

    def compute_cell_communication(
        self,
        mode: Literal["standard", "cell_type"] = "standard",
        model: str | None = None,
        n_permutations: int = 1000,
        seed: int = 42,
        test: Literal["parametric", "non-parametric", "both"] = "both",
        mean: Literal["algebraic", "geometric"] = "algebraic",
        layer_key_p_test: str | None = None,
        layer_key_np_test: str | None = None,
        device: str | None = None,
        verbose: bool = False,
    ) -> None:
        """Compute ligand-receptor / metabolite CCC scores.

        Requires ``compute_gene_pairs()`` to have been called.

        Parameters
        ----------
        mode
            ``"standard"`` for cell-type-agnostic analysis or ``"cell_type"``
            for cell-type-stratified analysis.
        model
            Expression normalization model.
        n_permutations
            Number of permutations for the non-parametric test.
        seed
            Random seed for reproducibility.
        test
            Statistical test to run.
        mean
            Averaging method for multi-gene interactions.
        layer_key_p_test
            Layer to use for the parametric test. Defaults to the layer set at
            construction time (``layer_key``).
        layer_key_np_test
            Layer to use for the non-parametric test. Defaults to the layer set at
            construction time (``layer_key``). Pass a different value from
            ``layer_key_p_test`` to use separate layers (e.g. raw counts vs log-norm).
        device
            PyTorch device string. ``None`` uses the underlying default.
        verbose
            Print progress.
        """
        self._require(STEP_CCC)
        if mode not in ("standard", "cell_type"):
            raise ValueError(f"mode must be 'standard' or 'cell_type', got '{mode}'")

        resolved_model = model if model is not None else "none"

        kwargs = {
            "adata": self._adata,
            "layer_key_p_test": (
                layer_key_p_test if layer_key_p_test is not None else self._layer_key
            ),
            "layer_key_np_test": (
                layer_key_np_test if layer_key_np_test is not None else self._layer_key
            ),
            "model": resolved_model,
            "M": n_permutations,
            "seed": seed,
            "test": test,
            "mean": mean,
            "device": device,
            "verbose": verbose,
        }
        if mode == "cell_type":
            _compute_ct_cell_communication(
                **kwargs,
                cell_type_key=self._adata.uns[HARREMAN_UNS_KEY].get("cell_type_key"),
            )
        else:
            _compute_cell_communication(**kwargs)

        self._ccc_mode = mode
        self._completed_steps.add(STEP_CCC)

    def compute_interacting_cell_scores(
        self,
        mode: Literal["standard", "cell_type"] = "standard",
        test: Literal["parametric", "non-parametric", "both"] = "both",
        restrict_significance: Literal["gene pairs", "metabolites", "both"] = "both",
        n_permutations: int = 1000,
        seed: int = 42,
        device: str | None = None,
        verbose: bool = False,
    ) -> None:
        """Compute per-cell interaction scores.

        Requires ``compute_cell_communication()`` to have been called.

        Parameters
        ----------
        mode
            ``"standard"`` or ``"cell_type"`` matching the CCC mode used.
        test
            Which test scores to compute.
        restrict_significance
            Which significant interactions to use from CCC results.
        n_permutations
            Permutations for non-parametric test.
        seed
            Random seed.
        device
            PyTorch device string.
        verbose
            Print progress.
        """
        self._require(STEP_ICS)
        if mode not in ("standard", "cell_type"):
            raise ValueError(f"mode must be 'standard' or 'cell_type', got '{mode}'")
        if self._ccc_mode is not None and mode != self._ccc_mode:
            import warnings

            warnings.warn(
                f"compute_interacting_cell_scores(mode='{mode}') does not match "
                f"compute_cell_communication(mode='{self._ccc_mode}'). "
                "Results may be inconsistent.",
                UserWarning,
                stacklevel=2,
            )

        if mode == "cell_type":
            _compute_ct_interacting_cell_scores(
                adata=self._adata,
                test=test,
                restrict_significance=restrict_significance,
                device=device,
                verbose=verbose,
            )
        else:
            _compute_interacting_cell_scores(
                adata=self._adata,
                test=test,
                restrict_significance=restrict_significance,
                M=n_permutations,
                seed=seed,
                device=device,
                verbose=verbose,
            )
        self._completed_steps.add(STEP_ICS)

    def select_significant_interactions(
        self,
        fdr_threshold: float = 0.05,
        test: Literal["parametric", "non-parametric"] = "parametric",
        use_fdr: bool = True,
        ct_aware: bool | None = None,
    ) -> None:
        """Filter to statistically significant interactions.

        Requires ``compute_cell_communication()`` to have been called.

        Parameters
        ----------
        fdr_threshold
            Significance cutoff.
        test
            Which test results to filter on.
        use_fdr
            Use FDR or raw p-values.
        ct_aware
            Whether to filter cell-type-aware results. Defaults to ``True`` when
            a cell type key was registered.
        """
        self._require(STEP_SIG)
        if ct_aware is None:
            ct_aware = self._ccc_mode == "cell_type"
        _select_significant_interactions(
            adata=self._adata,
            ct_aware=ct_aware,
            test=test,
            use_FDR=use_fdr,
            threshold=fdr_threshold,
        )
        self._completed_steps.add(STEP_SIG)


# ── Accessor classes ───────────────────────────────────────────────────────────


class _HarremanHsAccessor:
    """Hotspot analysis methods accessible as ``ha.hs.<method>``."""

    def __init__(self, ha: HarremanAnalysis) -> None:
        self._ha = ha

    def _resolve(self, adata: AnnData | None) -> AnnData:
        return adata if adata is not None else self._ha._adata

    def compute_local_autocorrelation(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.hotspot.local_autocorrelation import (
            compute_local_autocorrelation,
        )

        compute_local_autocorrelation(self._resolve(adata), **kwargs)

    def compute_local_correlation(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.hotspot.local_correlation import compute_local_correlation

        compute_local_correlation(self._resolve(adata), **kwargs)

    def create_modules(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.hotspot.modules import create_modules

        create_modules(self._resolve(adata), **kwargs)

    def calculate_module_scores(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.hotspot.modules import calculate_module_scores

        calculate_module_scores(self._resolve(adata), **kwargs)

    def calculate_super_module_scores(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.hotspot.modules import calculate_super_module_scores

        calculate_super_module_scores(self._resolve(adata), **kwargs)

    def compute_top_scoring_modules(self, adata: AnnData | None = None, **kwargs):
        from scvi.external.harreman.hotspot.modules import compute_top_scoring_modules

        return compute_top_scoring_modules(self._resolve(adata), **kwargs)

    def load_signatures(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.vision.signature import load_signatures

        load_signatures(self._resolve(adata), **kwargs)

    def compute_vision_signatures(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.vision.signature import compute_vision_signatures

        compute_vision_signatures(self._resolve(adata), **kwargs)

    def integrate_vision_hotspot_results(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.hotspot.modules import integrate_vision_hotspot_results

        integrate_vision_hotspot_results(self._resolve(adata), **kwargs)


class _HarremanTlAccessor:
    """Tool methods accessible as ``ha.tl.<method>``."""

    def __init__(self, ha: HarremanAnalysis) -> None:
        self._ha = ha

    def _resolve(self, adata: AnnData | None) -> AnnData:
        return adata if adata is not None else self._ha._adata

    def compute_knn_graph(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.tools.knn import compute_knn_graph

        compute_knn_graph(self._resolve(adata), **kwargs)

    def compute_interaction_module_correlation(
        self, adata: AnnData | None = None, **kwargs
    ) -> None:
        from scvi.external.harreman.tools.cell_communication import (
            compute_interaction_module_correlation,
        )

        compute_interaction_module_correlation(self._resolve(adata), **kwargs)

    def select_significant_interactions(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.tools.cell_communication import select_significant_interactions

        select_significant_interactions(self._resolve(adata), **kwargs)


class _HarremanPlAccessor:
    """Plotting methods accessible as ``ha.pl.<method>``."""

    def __init__(self, ha: HarremanAnalysis) -> None:
        self._ha = ha

    def _resolve(self, adata: AnnData | None) -> AnnData:
        return adata if adata is not None else self._ha._adata

    def local_correlation_plot(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.plots import local_correlation_plot

        local_correlation_plot(self._resolve(adata), **kwargs)

    def average_local_correlation_plot(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.plots import average_local_correlation_plot

        average_local_correlation_plot(self._resolve(adata), **kwargs)

    def module_score_correlation_plot(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.plots import module_score_correlation_plot

        module_score_correlation_plot(self._resolve(adata), **kwargs)

    def plot_interacting_cell_scores(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.plots import plot_interacting_cell_scores

        plot_interacting_cell_scores(self._resolve(adata), **kwargs)

    def plot_ct_interacting_cell_scores(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.plots import plot_ct_interacting_cell_scores

        plot_ct_interacting_cell_scores(self._resolve(adata), **kwargs)

    def plot_interaction_module_correlation(self, adata: AnnData | None = None, **kwargs) -> None:
        from scvi.external.harreman.plots import plot_interaction_module_correlation

        plot_interaction_module_correlation(self._resolve(adata), **kwargs)
