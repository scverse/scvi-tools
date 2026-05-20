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

| Model | Integration |
|-------|-------------|
| **DestVI** | `model.get_proportions(adata)` → cell-type layers |
| **SCVIVA** | `model.get_proportions()` + `model.get_latent_representation()` for KNN |
| **RESOLVI** | `model.get_denoised_expression()` as count layer |
| **None** | User provides pre-prepared adata (existing workflow) |

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
| `test_harreman_with_scviva` | SCVIVA | Cell-type proportions auto-extracted; `is_deconvolved=True` |
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
