# HarremanAnalysis OOP Design

**Date:** 2026-05-20
**Status:** Approved
**Branch:** add-harreman

---

## Decision Log

| # | Decision | Rationale |
|---|----------|-----------|
| D1 | Standalone analysis class (not mixin) | Harreman doesn't touch VAE internals. Works with any spatial adata. Model-agnostic by default. |
| D2 | AnnData required, model optional | Model is only needed for auto-extracting deconv proportions. Pre-prepared adata works without model. |
| D3 | New `_analysis.py` inside existing `harreman/` package | Co-located with tools it wraps. Existing files untouched. |
| D4 | Add `_constants.py` and `_results.py` alongside | Follows scvi-tools convention (resolvi, scviva each have `_constants.py`). |
| D5 | Extract model outputs at `__init__` time, drop model ref | Clean — no stale model references, no lazy extraction surprises. |
| D6 | Results in `adata.uns["harreman"]` + `HarremanResults` dataclass | Scanpy convention for persistence; dataclass for typed access. |
| D7 | Ordered step enforcement via `_completed_steps: set[str]` | Clear error messages when prerequisites are skipped. |
| D8 | No changes to existing harreman functional code | New class wraps existing functions. Enables behavior comparison. |

---

## Context

Harreman is a spatial metabolic cell-cell communication (CCC) toolkit that builds on Hotspot's local
autocorrelation statistics. It adds:
- Metabolite-level gene-pair integration via HarremanDB
- Cell-type-aware CCC inference
- Support for deconvolved spatial data (DestVI, SCVIVA)

The existing codebase is a flat collection of functional modules (`hotspot/`, `tools/`, `preprocessing/`).
This design wraps them into a production-ready OOP class compatible with scvi-tools conventions.

**Harreman vs Hotspot:** Hotspot is a building block (local autocorrelation, local correlation).
Harreman extends it with metabolite database integration and cell-type awareness.

---

## Target Models

| Model | What it provides | How HarremanAnalysis uses it |
|-------|-----------------|------------------------------|
| **DestVI** | Cell-type proportion matrix | `model.get_proportions()` → `(n_obs, n_labels)` → attach each cell-type as `adata.layers[ct]`, enables `mode="cell_type"` |
| **SCVIVA** | Cell-type labels in `adata.obs` + latent representation | `model.get_latent_representation()` for KNN coords; cell-type labels already in adata |
| **RESOLVI** | Denoised counts | `model.get_normalized_expression()` → attach as `adata.layers["denoised"]`, use as `layer_key` |
| **None** | — | User provides pre-prepared adata; existing pp.setup_anndata() workflow |

**Note:** SCVIVA does not have `get_proportions()`. It has `adata.obsm["niche_composition"]` (neighborhood composition) and cell-type labels in `adata.obs`. RESOLVI exposes `get_normalized_expression()` (not `get_denoised_expression()`) via `ResolVIPredictiveMixin`.

---

## File Structure

```
src/scvi/external/harreman/
  __init__.py          # updated: export HarremanAnalysis
  _analysis.py         # NEW: HarremanAnalysis class (~400 lines)
  _results.py          # NEW: HarremanResults dataclass
  _constants.py        # NEW: uns keys, registry keys, string constants
  hotspot/             # untouched
  tools/               # untouched
  preprocessing/       # untouched
  datasets/            # untouched

tests/external/
  test_harreman_analysis.py  # NEW: unit + integration tests
```

---

## Class API

```python
class HarremanAnalysis:
    """Downstream spatial metabolic cell-cell communication analysis.

    Parameters
    ----------
    adata
        Spatial AnnData. Must have neighbor weights in obsm.
    model
        Optional trained scvi spatial model (RESOLVI, SCVIVA, DestVI).
        DestVI: auto-calls get_proportions() and attaches as layers.
        RESOLVI/SCVIVA: uses latent coords for KNN if available.
        Model reference is dropped after init.
    layer_key
        Layer to use for counts. None = adata.X.

    Examples
    --------
    Without model:

    >>> ha = HarremanAnalysis(adata)
    >>> ha.setup(cell_type_key="cell_type", database_path="path/to/db.csv",
    ...          compute_neighbors_on_key="spatial")
    >>> ha.filter_genes()
    >>> ha.compute_gene_pairs()
    >>> ha.compute_cell_communication()
    >>> ha.select_significant_interactions()
    >>> results = ha.results

    With DestVI model:

    >>> ha = HarremanAnalysis(adata, model=destvi_model)
    >>> ha.setup(cell_type_key="cell_type", database_path="path/to/db.csv",
    ...          compute_neighbors_on_key="spatial")
    >>> ha.compute_cell_communication(mode="cell_type")
    """

    def __init__(
        self,
        adata: AnnData,
        model: BaseModelClass | None = None,
        layer_key: str | None = None,
    ) -> None: ...

    # ── Setup ─────────────────────────────────────────────────────────────
    def setup(
        self,
        cell_type_key: str,
        database_path: str | Path,
        compute_neighbors_on_key: str,
        spot_diameter: int = 10,
        sample_key: str | None = None,
        species: Literal["human", "mouse"] = "human",
    ) -> None:
        """Register adata with database and neighbor structure."""

    # ── Analysis steps ────────────────────────────────────────────────────
    def filter_genes(
        self,
        model: Literal["danb", "bernoulli", "normal", "none"] = "danb",
        feature_elimination: bool = False,
        threshold: float = 0.2,
        autocorrelation_filt: bool = False,
        expression_filt: bool = False,
        de_filt: bool = False,
        device: str | None = None,
    ) -> None:
        """Filter genes by sparsity, autocorrelation, and/or DE. Requires setup()."""

    def compute_gene_pairs(
        self,
        n_neighbors: int = 30,
        fdr_threshold: float = 0.05,
    ) -> None:
        """Compute local correlations for gene pairs in metabolite database. Requires setup()."""

    def compute_cell_communication(
        self,
        mode: Literal["standard", "cell_type"] = "standard",
        n_permutations: int = 200,
        fdr_threshold: float = 0.05,
    ) -> None:
        """Compute ligand-receptor CCC scores. Requires compute_gene_pairs()."""

    def compute_interacting_cell_scores(
        self,
        mode: Literal["standard", "cell_type"] = "standard",
    ) -> None:
        """Compute per-cell interaction scores. Requires compute_cell_communication()."""

    def select_significant_interactions(
        self,
        fdr_threshold: float = 0.05,
        min_cells: int = 5,
    ) -> None:
        """Filter to significant interactions. Requires compute_cell_communication()."""

    # ── Properties ────────────────────────────────────────────────────────
    @property
    def adata(self) -> AnnData: ...

    @property
    def results(self) -> HarremanResults:
        """Typed view over adata.uns['harreman']."""

    @property
    def is_deconvolved(self) -> bool:
        """True if model-extracted cell-type proportions are present."""

    @property
    def is_set_up(self) -> bool:
        """True if setup() has been called."""
```

---

## State Machine

Steps must run in order. `_completed_steps: set[str]` tracks progress.

```
[init] → setup() → filter_genes() → compute_gene_pairs()
                                          ↓
                              compute_cell_communication()
                                          ↓
                         compute_interacting_cell_scores()  (optional)
                                          ↓
                         select_significant_interactions()
```

`filter_genes()` is optional before `compute_gene_pairs()`.
`compute_interacting_cell_scores()` is optional.

Skipping a required prerequisite raises `RuntimeError`:
```
RuntimeError: compute_cell_communication() requires compute_gene_pairs() to be run first.
Call ha.compute_gene_pairs() before continuing.
```

---

## Results Storage

Results accumulate in `adata.uns["harreman"]` (survives `adata.write_h5ad()`):

```python
adata.uns["harreman"] = {
    "autocorrelation": pd.DataFrame,  # gene hotspot Z-scores
    "gene_pairs": pd.DataFrame,  # local correlations per metabolite pair
    "cell_communication": pd.DataFrame,  # LR CCC scores + p-values
    "ct_cell_communication": pd.DataFrame,  # cell-type-stratified CCC
    "interacting_cell_scores": ...,  # per-cell interaction scores
    "significant_interactions": pd.DataFrame,  # FDR-filtered
    "params": {  # reproducibility record
        "layer_key": ...,
        "model": ...,
        "species": ...,
        "filter_genes_kwargs": ...,
        "cell_communication_kwargs": ...,
    },
}
```

`ha.results` returns a `HarremanResults` dataclass that wraps these DataFrames (no copy).

---

## Constants (`_constants.py`)

```python
HARREMAN_UNS_KEY = "harreman"
HARREMAN_AUTOCORR_KEY = "autocorrelation"
HARREMAN_GENE_PAIRS_KEY = "gene_pairs"
HARREMAN_CCC_KEY = "cell_communication"
HARREMAN_CT_CCC_KEY = "ct_cell_communication"
HARREMAN_ICS_KEY = "interacting_cell_scores"
HARREMAN_SIG_KEY = "significant_interactions"
HARREMAN_PARAMS_KEY = "params"

STEP_SETUP = "setup"
STEP_FILTER = "filter_genes"
STEP_GENE_PAIRS = "compute_gene_pairs"
STEP_CCC = "compute_cell_communication"
STEP_ICS = "compute_interacting_cell_scores"
STEP_SIG = "select_significant_interactions"

SUPPORTED_MODELS = ("RESOLVI", "SCVIVA", "DestVI")
```

---

## Results Dataclass (`_results.py`)

```python
@dataclass
class HarremanResults:
    """Typed view over adata.uns['harreman']."""

    autocorrelation: pd.DataFrame | None
    gene_pairs: pd.DataFrame | None
    cell_communication: pd.DataFrame | None
    ct_cell_communication: pd.DataFrame | None
    interacting_cell_scores: Any | None
    significant_interactions: pd.DataFrame | None
    params: dict

    @classmethod
    def from_uns(cls, uns: dict) -> "HarremanResults":
        """Construct from adata.uns['harreman']."""
```

---

## Testing Plan

### Unit tests (fast, no model training)
File: `tests/external/test_harreman_analysis.py`

| Test | Coverage |
|------|----------|
| `test_init_no_model` | adata-only init, adata stored correctly |
| `test_init_with_destvi_mock` | proportions extracted, layers attached (MagicMock) |
| `test_init_with_resolvi_mock` | denoised counts extracted (MagicMock) |
| `test_init_unsupported_model` | ValueError with clear message |
| `test_setup` | database loaded, obs/uns fields set |
| `test_filter_genes` | prerequisite check, gene mask applied |
| `test_compute_gene_pairs` | output shapes, uns keys populated |
| `test_compute_cell_communication_standard` | results DataFrame non-empty |
| `test_compute_cell_communication_cell_type` | ct-mode results |
| `test_compute_interacting_cell_scores` | per-cell scores |
| `test_select_significant_interactions` | FDR filtering applied |
| `test_prerequisite_order_enforced` | RuntimeError if steps skipped |
| `test_results_property` | HarremanResults wraps adata.uns correctly |
| `test_is_deconvolved_property` | True when proportions present |

**Mock strategy:** Synthetic AnnData (20 cells, 50 genes, 2 cell types). Mock model via
`unittest.mock.MagicMock` with `get_proportions()` / `get_denoised_expression()` returning
fixture DataFrames. No real spatial data required for unit tests.

### Integration tests (real model training, `@pytest.mark.optional`)
File: `tests/external/test_harreman_analysis.py` (same file, `@pytest.mark.optional`)

All three spatial models are tested end-to-end: train minimal model (1-2 epochs on
`synthetic_iid(generate_coordinates=True)`), pass to `HarremanAnalysis`, run full pipeline.

| Test | Model | Key assertion |
|------|-------|---------------|
| `test_harreman_with_resolvi` | RESOLVI | Denoised counts attached as layer; full pipeline runs |
| `test_harreman_with_scviva` | SCVIVA | Latent representation used for KNN; cell-type labels from `adata.obs` |
| `test_harreman_with_destvi` | DestVI | `get_proportions()` called once; cell-type layers match shape |

**Fixture reuse:** Each integration test uses a session-scoped fixture that trains the model
once per test session (following the pattern in `tests/external/resolvi/test_resolvi.py`
and `tests/model/test_destvi.py`).

```python
@pytest.fixture(scope="session")
def resolvi_model():
    adata = synthetic_iid(generate_coordinates=True, n_regions=5)
    adata.obsm["X_spatial"] = adata.obsm["coordinates"]
    RESOLVI.setup_anndata(adata)
    model = RESOLVI(adata)
    model.train(max_epochs=2)
    return model, adata


@pytest.fixture(scope="session")
def destvi_model():
    # CondSCVI → DestVI as in tests/model/test_destvi.py
    dataset = synthetic_iid(n_labels=5)
    sc_model = CondSCVI(dataset, n_latent=2)
    sc_model.train(1)
    DestVI.setup_anndata(dataset)
    spatial_model = DestVI.from_rna_model(dataset, sc_model)
    spatial_model.train(max_epochs=1)
    return spatial_model, dataset


@pytest.fixture(scope="session")
def scviva_model():
    adata = synthetic_iid(batch_size=256, n_genes=100, generate_coordinates=True)
    SCVIVA.preprocessing_anndata(adata, k_nn=5, ...)
    SCVIVA.setup_anndata(adata, layer="counts", ...)
    model = SCVIVA(adata)
    model.train(max_epochs=2)
    return model, adata
```

---

## Open Questions (for future)

- CellAssign support: add as a 4th model type if cell type annotations needed from it
- Export / plotting: `ha.plot_interactions()` could be a future addition
- Multi-sample support: `sample_key` is accepted but multi-sample logic TBD

---

## Future Work: visionpy Roadmap

### Background

[visionpy](https://github.com/YosefLab/visionpy) is Adam's Python reimplementation of the
[R VISION package](https://github.com/YosefLab/VISION). Several Harreman functions (gene
signature scoring, local autocorrelation via Geary-C) overlap with VISION's analysis surface.
visionpy is explicitly "NOT YET READY FOR USE" and is missing significant R functionality.

### Current State (as of 2026-05-20)

**visionpy has:**
- `signature.py` — gene signature scoring
- `diffexp.py` — differential expression
- `anndata.py` — AnnData I/O helpers
- `api.py` + `blueprint.py` + Flask web layer — interactive report server

**R VISION has (not yet in visionpy):**

| R file | Functionality |
|--------|--------------|
| `Microclusters.R` | Microcluster pooling — group similar cells before analysis to handle large datasets |
| `FitchParsimony.R` + `TreeBasedMethods.R` | **PhyloVISION** — tree-based signature analysis using Fitch parsimony over cell lineage trees |
| `Projections.R` | Unified dimensionality reduction framework (PCA, UMAP, t-SNE, trajectory coords) |
| `methods-Trajectory.R` + `methods-TrajectoryProjection.R` | Trajectory integration — score signatures along pseudotime axes |
| `NormalizationMethods.R` | Library-size normalization, log-transform, FNR correction |
| `RcppExports.R` | C++ speed layer (Geary-C statistic, permutation tests) |
| `Filters.R` | Gene/cell filtering prior to analysis |

### Architectural Decision: restructure visionpy for scvi-tools compatibility

**Recommendation: yes, restructure.** Rationale:

1. visionpy is pre-1.0 with no stable API — now is the cheapest time to restructure.
2. VISION's natural scvi-tools counterpart is a downstream analysis class, analogous to
   `HarremanAnalysis` — it takes a trained model + AnnData and wraps the analysis pipeline.
3. scvi-tools already provides the spatial/embedding models (SCVIVA, RESOLVI, DestVI) that
   VISION would annotate. A `VisionAnalysis` class could accept these models directly.
4. Hotspot/Harreman is already being restructured this way — visionpy should follow the same
   pattern for consistency.

**Proposed class:** `VisionAnalysis(adata, model=None, signatures=None)`

### Planned additions (priority order)

#### 1. Hotspot / Harreman integration (fast win)

- `VisionAnalysis` consumes `HarremanAnalysis.results` as an optional input.
- Autocorrelation results from Hotspot (`gene_autocorrelation_results`) surface as signature
  consistency scores within the VISION report.
- Shared KNN graph: build once in `HarremanAnalysis.setup()`, reuse in `VisionAnalysis`.
- Signature-level Geary-C scores replace per-signature permutation tests where Hotspot
  Z-scores are available.

#### 2. Microcluster pooling

- Port `Microclusters.R` logic to Python/NumPy.
- Expose as `VisionAnalysis.pool_microclusters(n_cells=10)` — optional preprocessing step
  before signature scoring to enable scalability to millions of cells.
- Store pooled representation in `adata.uns["vision"]["microclusters"]`.
- All downstream analysis runs on pooled cells; results are unpooled back to cell level.

#### 3. PhyloVISION integration

- Port `FitchParsimony.R` + `TreeBasedMethods.R`.
- Input: a Newick/ETE3 lineage tree over cells (e.g., from scLineage or cassiopeia).
- `VisionAnalysis.compute_phylo_scores(tree)` — run Fitch parsimony on each signature over
  the tree to score lineage-consistency.
- Store in `adata.uns["vision"]["phylo_scores"]`.
- Optional step, gated by whether a tree is provided.

#### 4. Additional parity items

- `Projections.R` → `VisionAnalysis.add_projection(key, coords)` — register arbitrary
  low-dimensional embeddings to score signatures against.
- Trajectory support → `VisionAnalysis.add_trajectory(pseudotime_key)` — score signatures
  along pseudotime.
- `NormalizationMethods.R` → normalization handled upstream by scvi-tools models; only need
  a thin wrapper for users supplying raw counts without a model.
- `Filters.R` → delegate to `HarremanAnalysis.filter_genes()` or scanpy's `sc.pp.filter_genes`.

#### 5. PyTorch scalability

- Replace C++ Rcpp layer (Geary-C, permutation tests) with PyTorch batched operations.
- Target: handle 1M+ cells via GPU-accelerated sparse matrix ops.
- Microcluster pooling reduces effective N before the bottleneck computations.
- KNN graph construction via `faiss` or `pynndescent` (already used in scvi-tools).
- Permutation tests via `torch.randperm` on GPU, parallelized across signatures.

#### 6. scvi-tools integration pattern

```
VisionAnalysis(
    adata,
    model=scviva_model,        # optional — provides latent coords + cell-type labels
    signatures=sig_library,    # MSigDB-format dict or .gmt path
    harreman=ha,               # optional HarremanAnalysis instance to reuse KNN + autocorr
)
```

- Same model dispatch table as `HarremanAnalysis`: SCVIVA → latent coords, RESOLVI → denoised
  expression, DestVI → proportions.
- Results in `adata.uns["vision"]` following same pattern as `adata.uns["harreman"]`.
- `VisionResults` dataclass mirrors `HarremanResults`.

### Proposed file structure

```
src/scvi/external/vision/        # new package, parallel to harreman/
  __init__.py
  _analysis.py                   # VisionAnalysis class
  _results.py                    # VisionResults dataclass
  _constants.py                  # uns keys, step names
  _microclusters.py              # microcluster pooling
  _phylo.py                      # PhyloVISION (Fitch parsimony)
  _projections.py                # projection registry
  _scoring.py                    # signature scoring (port of visionpy signature.py)
  _server.py                     # Flask web report (port of visionpy blueprint.py)
```

visionpy upstream repo: keep or archive. If restructured into scvi-tools, the standalone
visionpy package becomes redundant; coordinate with Adam before deprecating.

### Dependencies not yet in scvi-tools

| Dependency | Purpose | Notes |
|------------|---------|-------|
| `ete3` or `cassiopeia` | Lineage tree parsing for PhyloVISION | optional extra |
| `flask` | Interactive web report | already in visionpy |
| `faiss-cpu` / `faiss-gpu` | Fast KNN for large datasets | check if already present |
