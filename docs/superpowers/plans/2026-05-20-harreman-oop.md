# HarremanAnalysis OOP Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wrap the existing Harreman functional tools into a production-ready `HarremanAnalysis` standalone class compatible with scvi-tools conventions, with integration for RESOLVI, SCVIVA, and DestVI.

**Architecture:** `HarremanAnalysis` is a stateful analysis class (like `PosteriorPredictiveCheck`) that wraps the existing functional API in `external/harreman/tools/`, `hotspot/`, and `preprocessing/`. It stores results in `adata.uns["harreman"]` following scanpy convention and validates step ordering via `_completed_steps`. Model integration (DestVI/RESOLVI/SCVIVA) extracts needed outputs at `__init__` time and attaches them as layers or obsm entries; the model reference is then dropped.

**Tech Stack:** Python 3.12+, AnnData, pandas, numpy, torch, pytest, unittest.mock. All underlying compute reuses the existing harreman functional code untouched. Use `X | Y` unions, lowercase generics (`list[str]`), keep `from __future__ import annotations` per repo convention.

**Spec:** `docs/superpowers/specs/2026-05-20-harreman-oop-design.md`

**Note on database:** HarremanDB is downloaded automatically from S3 via `pooch` inside `preprocessing.database.extract_interaction_db`. There is no `database_path` parameter — callers pass `species` and `database` (type). Unit tests mock this function to avoid network calls.

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `src/scvi/external/harreman/_constants.py` | Create | All string constants: uns keys, step names, model type strings |
| `src/scvi/external/harreman/_results.py` | Create | `HarremanResults` dataclass wrapping `adata.uns["harreman"]` |
| `src/scvi/external/harreman/_analysis.py` | Create | `HarremanAnalysis` class (~450 lines) |
| `src/scvi/external/harreman/__init__.py` | Modify | Export `HarremanAnalysis` |
| `tests/external/harreman/test_harreman_analysis.py` | Create | Unit tests + `@pytest.mark.optional` integration tests |

**Do not modify** any file in `hotspot/`, `tools/`, `preprocessing/`, or `datasets/`.

---

## Task 1: `_constants.py`

**Files:**
- Create: `src/scvi/external/harreman/_constants.py`

- [ ] **Step 1: Create the constants file**

```python
# src/scvi/external/harreman/_constants.py
from __future__ import annotations

# ── adata.uns keys ────────────────────────────────────────────────────────────
HARREMAN_UNS_KEY = "harreman"
HARREMAN_AUTOCORR_KEY = "autocorrelation"
HARREMAN_GENE_PAIRS_KEY = "gene_pairs_results"
HARREMAN_CCC_KEY = "cell_communication"
HARREMAN_CT_CCC_KEY = "ct_cell_communication"
HARREMAN_ICS_KEY = "interacting_cell_scores"
HARREMAN_SIG_KEY = "significant_interactions"
HARREMAN_PARAMS_KEY = "params"

# ── Step names (used in _completed_steps set) ─────────────────────────────────
STEP_SETUP = "setup"
STEP_FILTER = "filter_genes"
STEP_GENE_PAIRS = "compute_gene_pairs"
STEP_CCC = "compute_cell_communication"
STEP_ICS = "compute_interacting_cell_scores"
STEP_SIG = "select_significant_interactions"

# ── Model type strings ────────────────────────────────────────────────────────
MODEL_DESTVI = "DestVI"
MODEL_RESOLVI = "RESOLVI"
MODEL_SCVIVA = "SCVIVA"

SUPPORTED_MODELS = (MODEL_DESTVI, MODEL_RESOLVI, MODEL_SCVIVA)

# ── Internal adata keys written during model integration ─────────────────────
HARREMAN_DENOISED_LAYER = "harreman_denoised"
HARREMAN_LATENT_OBSM = "X_harreman_latent"
```

- [ ] **Step 2: Verify no import errors**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -c "from scvi.external.harreman._constants import HARREMAN_UNS_KEY, SUPPORTED_MODELS; print('OK')"
```

Expected: `OK`

---

## Task 2: `_results.py` + tests

**Files:**
- Create: `src/scvi/external/harreman/_results.py`
- Create: `tests/external/harreman/test_harreman_analysis.py` (first section)

- [ ] **Step 1: Write failing tests for HarremanResults**

Create `tests/external/harreman/test_harreman_analysis.py`:

```python
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from unittest.mock import MagicMock, patch

from scvi.external.harreman._constants import (
    HARREMAN_AUTOCORR_KEY,
    HARREMAN_CCC_KEY,
    HARREMAN_CT_CCC_KEY,
    HARREMAN_GENE_PAIRS_KEY,
    HARREMAN_ICS_KEY,
    HARREMAN_PARAMS_KEY,
    HARREMAN_SIG_KEY,
    HARREMAN_UNS_KEY,
    STEP_CCC,
    STEP_FILTER,
    STEP_GENE_PAIRS,
    STEP_SETUP,
)
from scvi.external.harreman._results import HarremanResults

# ── Shared fixtures ────────────────────────────────────────────────────────────


@pytest.fixture
def adata_spatial():
    n_obs, n_vars = 50, 20
    rng = np.random.default_rng(42)
    X = rng.poisson(1.0, size=(n_obs, n_vars)).astype(float)
    obs = pd.DataFrame(
        {"cell_type": pd.Categorical(np.tile(["TypeA", "TypeB"], n_obs // 2))},
        index=[f"cell{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene{i}" for i in range(n_vars)])
    obsm = {"spatial": rng.random((n_obs, 2)) * 100}
    return AnnData(X=X, obs=obs, var=var, obsm=obsm)


@pytest.fixture
def uns_with_results():
    df = pd.DataFrame({"score": [0.1, 0.2]})
    return {
        HARREMAN_AUTOCORR_KEY: df.copy(),
        HARREMAN_GENE_PAIRS_KEY: df.copy(),
        HARREMAN_CCC_KEY: df.copy(),
        HARREMAN_CT_CCC_KEY: None,
        HARREMAN_ICS_KEY: None,
        HARREMAN_SIG_KEY: df.copy(),
        HARREMAN_PARAMS_KEY: {"species": "human"},
    }


# ── HarremanResults tests ──────────────────────────────────────────────────────


def test_results_from_uns(uns_with_results):
    r = HarremanResults.from_uns(uns_with_results)
    assert isinstance(r.autocorrelation, pd.DataFrame)
    assert isinstance(r.cell_communication, pd.DataFrame)
    assert r.ct_cell_communication is None
    assert r.params == {"species": "human"}


def test_results_from_uns_missing_key_raises(uns_with_results):
    del uns_with_results[HARREMAN_PARAMS_KEY]
    with pytest.raises(KeyError):
        HarremanResults.from_uns(uns_with_results)


def test_results_repr_does_not_raise(uns_with_results):
    r = HarremanResults.from_uns(uns_with_results)
    repr(r)  # should not raise
```

- [ ] **Step 2: Run tests — expect ImportError (class doesn't exist)**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py::test_results_from_uns -v 2>&1 | tail -10
```

Expected: `ImportError` or `ModuleNotFoundError` for `_results`

- [ ] **Step 3: Implement `_results.py`**

```python
# src/scvi/external/harreman/_results.py
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from scvi.external.harreman._constants import (
    HARREMAN_AUTOCORR_KEY,
    HARREMAN_CCC_KEY,
    HARREMAN_CT_CCC_KEY,
    HARREMAN_GENE_PAIRS_KEY,
    HARREMAN_ICS_KEY,
    HARREMAN_PARAMS_KEY,
    HARREMAN_SIG_KEY,
)

if TYPE_CHECKING:
    import pandas as pd


@dataclass
class HarremanResults:
    """Typed view over ``adata.uns['harreman']``.

    All DataFrame fields are direct references to the stored objects — no copy.
    """

    autocorrelation: pd.DataFrame | None
    gene_pairs: pd.DataFrame | None
    cell_communication: pd.DataFrame | None
    ct_cell_communication: pd.DataFrame | None
    interacting_cell_scores: Any | None
    significant_interactions: pd.DataFrame | None
    params: dict

    @classmethod
    def from_uns(cls, uns: dict) -> HarremanResults:
        """Construct from ``adata.uns['harreman']`` dict."""
        return cls(
            autocorrelation=uns.get(HARREMAN_AUTOCORR_KEY),
            gene_pairs=uns.get(HARREMAN_GENE_PAIRS_KEY),
            cell_communication=uns.get(HARREMAN_CCC_KEY),
            ct_cell_communication=uns.get(HARREMAN_CT_CCC_KEY),
            interacting_cell_scores=uns.get(HARREMAN_ICS_KEY),
            significant_interactions=uns.get(HARREMAN_SIG_KEY),
            params=uns[HARREMAN_PARAMS_KEY],
        )
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -k "results" -v 2>&1 | tail -15
```

Expected: 3 PASSED

---

## Task 3: `HarremanAnalysis` skeleton — init, properties, prerequisite checking

**Files:**
- Create: `src/scvi/external/harreman/_analysis.py`
- Modify: `tests/external/harreman/test_harreman_analysis.py` (add init/property tests)

- [ ] **Step 1: Add failing tests for init and properties**

Append to `tests/external/harreman/test_harreman_analysis.py`:

```python
from scvi.external.harreman._analysis import HarremanAnalysis

# ── Init / basic validation ────────────────────────────────────────────────────


def test_init_stores_adata(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    assert ha.adata is adata_spatial


def test_init_non_anndata_raises():
    with pytest.raises(TypeError, match="AnnData"):
        HarremanAnalysis("not_an_adata")


def test_init_is_not_deconvolved_by_default(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    assert ha.is_deconvolved is False


def test_init_is_not_set_up_by_default(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    assert ha.is_set_up is False


def test_results_before_any_steps_raises(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    with pytest.raises(RuntimeError, match="setup"):
        _ = ha.results


def test_filter_genes_before_setup_raises(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    with pytest.raises(RuntimeError, match="setup"):
        ha.filter_genes()


def test_compute_gene_pairs_before_setup_raises(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    with pytest.raises(RuntimeError, match="setup"):
        ha.compute_gene_pairs()


def test_compute_ccc_before_gene_pairs_raises(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    ha._completed_steps.add(STEP_SETUP)  # bypass setup check
    with pytest.raises(RuntimeError, match="compute_gene_pairs"):
        ha.compute_cell_communication()


def test_unsupported_model_raises(adata_spatial):
    bad_model = MagicMock()
    bad_model.__class__.__name__ = "NotAKnownModel"
    with pytest.raises(ValueError, match="Unsupported model"):
        HarremanAnalysis(adata_spatial, model=bad_model)
```

- [ ] **Step 2: Run tests — expect ImportError**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -k "init or raises or set_up or deconvolved" -v 2>&1 | tail -15
```

Expected: `ImportError` for `_analysis`

- [ ] **Step 3: Implement `_analysis.py` skeleton**

```python
# src/scvi/external/harreman/_analysis.py
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Literal

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

if TYPE_CHECKING:
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

        if model is not None:
            self._apply_model_integration(model, adata)

        self._adata.uns.setdefault(HARREMAN_UNS_KEY, {HARREMAN_PARAMS_KEY: {}})

    # ── Model integration ──────────────────────────────────────────────────────

    def _apply_model_integration(self, model: BaseModelClass, adata: AnnData) -> None:
        model_type = type(model).__name__
        if model_type not in SUPPORTED_MODELS:
            raise ValueError(
                f"Unsupported model type '{model_type}'. "
                f"Supported: {SUPPORTED_MODELS}"
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
            "HarremanAnalysis: attached %d DestVI cell-type layers.",
            len(proportions.columns),
        )

    @staticmethod
    def _extract_resolvi_outputs(model: BaseModelClass, adata: AnnData) -> None:
        """Attach denoised counts from RESOLVI as a layer."""
        denoised = model.get_normalized_expression()  # DataFrame (n_obs, n_vars)
        import numpy as np
        import pandas as pd

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

    # ── Prerequisite checking ──────────────────────────────────────────────────

    def _require(self, step: str) -> None:
        prereq = _STEP_PREREQUISITES.get(step)
        if prereq is not None and prereq not in self._completed_steps:
            prereq_method = prereq.replace("_", " ")
            step_method = step.replace("_", " ")
            raise RuntimeError(
                f"{step_method}() requires {prereq_method}() to be run first. "
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
            raise RuntimeError(
                "No results available. Run ha.setup() first, then analysis steps."
            )
        return HarremanResults.from_uns(self._adata.uns[HARREMAN_UNS_KEY])

    @property
    def is_deconvolved(self) -> bool:
        """True if model-extracted cell-type proportion layers are present."""
        return self._is_deconvolved

    @property
    def is_set_up(self) -> bool:
        """True if setup() has been called."""
        return STEP_SETUP in self._completed_steps
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -k "init or raises or set_up or deconvolved" -v 2>&1 | tail -20
```

Expected: all PASSED (ImportError gone, RuntimeError/ValueError/TypeError checks pass)

---

## Task 4: `setup()` method + tests

**Files:**
- Modify: `src/scvi/external/harreman/_analysis.py` (add `setup()`)
- Modify: `tests/external/harreman/test_harreman_analysis.py` (add setup tests)

- [ ] **Step 1: Add failing tests for setup()**

Append to test file:

```python
# ── setup() tests ─────────────────────────────────────────────────────────────


def _make_mock_extract_db(adata_target):
    """Return a patch that fakes extract_interaction_db side-effects."""
    import numpy as np

    def _fake_extract(adata, **kwargs):
        n_vars = adata.n_vars
        import pandas as pd

        db_df = pd.DataFrame(
            np.zeros((n_vars, 3)),
            index=adata.var_names,
            columns=["metab_A", "metab_B", "metab_C"],
        )
        adata.varm["database"] = db_df
        adata.uns["database_varm_key"] = "database"
        adata.uns["database"] = "transporter"
        adata.uns["metabolite_database"] = pd.DataFrame()
        adata.uns["num_metabolites"] = 3
        adata.uns["importer"] = pd.DataFrame()
        adata.uns["exporter"] = pd.DataFrame()
        adata.uns["import_export"] = pd.DataFrame()

    return _fake_extract


def test_setup_sets_is_set_up(adata_spatial, monkeypatch):
    import scvi.external.harreman._analysis as _mod

    monkeypatch.setattr(
        _mod, "_extract_interaction_db", _make_mock_extract_db(adata_spatial)
    )
    ha = HarremanAnalysis(adata_spatial)
    ha.setup(compute_neighbors_on_key="spatial", n_neighbors=5)
    assert ha.is_set_up is True


def test_setup_builds_knn_graph(adata_spatial, monkeypatch):
    import scvi.external.harreman._analysis as _mod

    monkeypatch.setattr(
        _mod, "_extract_interaction_db", _make_mock_extract_db(adata_spatial)
    )
    ha = HarremanAnalysis(adata_spatial)
    ha.setup(compute_neighbors_on_key="spatial", n_neighbors=5)
    assert "weights" in adata_spatial.obsp


def test_setup_writes_params(adata_spatial, monkeypatch):
    import scvi.external.harreman._analysis as _mod

    monkeypatch.setattr(
        _mod, "_extract_interaction_db", _make_mock_extract_db(adata_spatial)
    )
    ha = HarremanAnalysis(adata_spatial)
    ha.setup(compute_neighbors_on_key="spatial", n_neighbors=5, species="mouse")
    params = adata_spatial.uns[HARREMAN_UNS_KEY][HARREMAN_PARAMS_KEY]
    assert params["species"] == "mouse"


def test_setup_second_call_overwrites(adata_spatial, monkeypatch):
    import scvi.external.harreman._analysis as _mod

    monkeypatch.setattr(
        _mod, "_extract_interaction_db", _make_mock_extract_db(adata_spatial)
    )
    ha = HarremanAnalysis(adata_spatial)
    ha.setup(compute_neighbors_on_key="spatial", n_neighbors=5, species="human")
    ha.setup(compute_neighbors_on_key="spatial", n_neighbors=5, species="mouse")
    assert (
        adata_spatial.uns[HARREMAN_UNS_KEY][HARREMAN_PARAMS_KEY]["species"] == "mouse"
    )
```

- [ ] **Step 2: Run tests — expect FAIL (setup not implemented)**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -k "setup" -v 2>&1 | tail -15
```

Expected: AttributeError — `setup` not found

- [ ] **Step 3: Add `setup()` to `_analysis.py`**

Add to imports at top of `_analysis.py`:

```python
from scvi.external.harreman.preprocessing.database import (
    extract_interaction_db as _extract_interaction_db,
)
from scvi.external.harreman.tools.knn import compute_knn_graph as _compute_knn_graph
```

Add the method inside `HarremanAnalysis` (after `_extract_scviva_outputs`):

```python
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
) -> None:
    """Register adata with the HarremanDB and build the spatial graph.

    Parameters
    ----------
    compute_neighbors_on_key
        Key in ``adata.obsm`` to use for KNN computation (e.g. ``"spatial"``).
        For SCVIVA, pass ``HARREMAN_LATENT_OBSM``.
    species
        ``"human"`` or ``"mouse"``.
    database
        Which database to load: ``"transporter"``, ``"LR"``, or ``"both"``.
    n_neighbors
        Number of nearest neighbors for the spatial graph. Either this or
        ``neighborhood_radius`` must be provided.
    neighborhood_radius
        Radius for neighborhood graph. Alternative to ``n_neighbors``.
    cell_type_key
        Key in ``adata.obs`` for cell types. Required for ``mode="cell_type"``
        in downstream steps.
    sample_key
        Key in ``adata.obs`` for sample/batch. Used in multi-sample analysis.
    spot_diameter
        Spot diameter of the spatial technology (e.g. 10 for Visium).
    extracellular_only
        Restrict HarremanDB to extracellular metabolites only.
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
    )
    if cell_type_key is not None:
        self._adata.uns["cell_type_key"] = cell_type_key
    if sample_key is not None:
        self._adata.uns["sample_key"] = sample_key
    if spot_diameter is not None:
        self._adata.uns["spot_diameter"] = spot_diameter

    self._adata.uns[HARREMAN_UNS_KEY][HARREMAN_PARAMS_KEY].update(
        {
            "species": species,
            "database": database,
            "layer_key": self._layer_key,
            "compute_neighbors_on_key": compute_neighbors_on_key,
            "cell_type_key": cell_type_key,
        }
    )
    self._completed_steps.add(STEP_SETUP)
    logger.info(
        "HarremanAnalysis: setup complete (species=%s, database=%s).", species, database
    )
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -k "setup" -v 2>&1 | tail -15
```

Expected: 4 PASSED

---

## Task 5: `filter_genes()` + tests

**Files:**
- Modify: `src/scvi/external/harreman/_analysis.py`
- Modify: `tests/external/harreman/test_harreman_analysis.py`

- [ ] **Step 1: Add failing test**

Append to test file:

```python
# ── filter_genes() tests ───────────────────────────────────────────────────────


def _setup_ha_no_network(adata, monkeypatch, cell_type_key=None):
    """Helper: create HarremanAnalysis with setup() mocked to avoid S3."""
    import scvi.external.harreman._analysis as _mod

    monkeypatch.setattr(_mod, "_extract_interaction_db", _make_mock_extract_db(adata))
    ha = HarremanAnalysis(adata)
    ha.setup(
        compute_neighbors_on_key="spatial",
        n_neighbors=5,
        cell_type_key=cell_type_key,
    )
    return ha


def test_filter_genes_marks_step_complete(adata_spatial, monkeypatch):
    ha = _setup_ha_no_network(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._apply_gene_filtering") as mock_filt:
        mock_filt.return_value = None
        ha.filter_genes()
    assert STEP_FILTER in ha._completed_steps


def test_filter_genes_calls_underlying_function(adata_spatial, monkeypatch):
    ha = _setup_ha_no_network(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._apply_gene_filtering") as mock_filt:
        mock_filt.return_value = None
        ha.filter_genes(feature_elimination=True, threshold=0.3)
    mock_filt.assert_called_once()
    _, kwargs = mock_filt.call_args
    assert kwargs.get("feature_elimination") is True
    assert kwargs.get("threshold") == 0.3
```

- [ ] **Step 2: Run tests — expect FAIL**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -k "filter_genes" -v 2>&1 | tail -10
```

Expected: AttributeError — `filter_genes` not found

- [ ] **Step 3: Add `filter_genes()` to `_analysis.py`**

Add to imports:

```python
from scvi.external.harreman.tools.cell_communication import (
    apply_gene_filtering as _apply_gene_filtering,
)
```

Add method inside class:

```python
# ── Analysis steps ─────────────────────────────────────────────────────────


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
    """Filter genes before gene pair computation.

    Parameters
    ----------
    model
        Expression model for autocorrelation-based filtering.
    feature_elimination
        Filter genes expressed in fewer than ``threshold`` fraction of cells.
    threshold
        Minimum fraction of cells expressing a gene (used with ``feature_elimination``).
    autocorrelation_filt
        Keep only genes with significant local autocorrelation.
    expression_filt
        Keep only genes expressed in each cell type.
    de_filt
        Keep only genes differentially expressed between cell types.
    device
        PyTorch device string for GPU computation. ``None`` uses auto-selection.
    """
    self._require(STEP_FILTER)
    _apply_gene_filtering(
        adata=self._adata,
        layer_key=self._layer_key,
        cell_type_key=self._adata.uns.get("cell_type_key"),
        model=model,
        feature_elimination=feature_elimination,
        threshold=threshold,
        autocorrelation_filt=autocorrelation_filt,
        expression_filt=expression_filt,
        de_filt=de_filt,
        umi_counts_obs_key=None,
        device=device,
    )
    self._completed_steps.add(STEP_FILTER)
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -k "filter_genes" -v 2>&1 | tail -10
```

Expected: 2 PASSED

---

## Task 6: `compute_gene_pairs()` + tests

**Files:**
- Modify: `src/scvi/external/harreman/_analysis.py`
- Modify: `tests/external/harreman/test_harreman_analysis.py`

- [ ] **Step 1: Add failing test**

Append to test file:

```python
# ── compute_gene_pairs() tests ─────────────────────────────────────────────────


def test_compute_gene_pairs_marks_step_complete(adata_spatial, monkeypatch):
    ha = _setup_ha_no_network(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._compute_gene_pairs") as mock_gp:
        mock_gp.return_value = None
        ha.compute_gene_pairs()
    assert STEP_GENE_PAIRS in ha._completed_steps


def test_compute_gene_pairs_passes_layer_key(adata_spatial, monkeypatch):
    import scvi.external.harreman._analysis as _mod

    monkeypatch.setattr(
        _mod, "_extract_interaction_db", _make_mock_extract_db(adata_spatial)
    )
    ha = HarremanAnalysis(adata_spatial, layer_key="my_layer")
    ha.setup(compute_neighbors_on_key="spatial", n_neighbors=5)
    with patch("scvi.external.harreman._analysis._compute_gene_pairs") as mock_gp:
        mock_gp.return_value = None
        ha.compute_gene_pairs()
    _, kwargs = mock_gp.call_args
    assert kwargs.get("layer_key") == "my_layer"
```

- [ ] **Step 2: Run tests — expect FAIL**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -k "gene_pairs" -v 2>&1 | tail -10
```

Expected: AttributeError

- [ ] **Step 3: Add `compute_gene_pairs()` to `_analysis.py`**

Add to imports:

```python
from scvi.external.harreman.tools.cell_communication import (
    compute_gene_pairs as _compute_gene_pairs,
)
```

Add method:

```python
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
        Fix one cell type across all pairs (e.g. restrict to sender).
    verbose
        Print progress messages.
    """
    self._require(STEP_GENE_PAIRS)
    _compute_gene_pairs(
        adata=self._adata,
        layer_key=self._layer_key,
        cell_type_key=self._adata.uns.get("cell_type_key"),
        cell_type_pairs=cell_type_pairs,
        ct_specific=ct_specific,
        fix_ct=fix_ct,
        verbose=verbose,
    )
    self._completed_steps.add(STEP_GENE_PAIRS)
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -k "gene_pairs" -v 2>&1 | tail -10
```

Expected: 2 PASSED

---

## Task 7: `compute_cell_communication()` + tests

**Files:**
- Modify: `src/scvi/external/harreman/_analysis.py`
- Modify: `tests/external/harreman/test_harreman_analysis.py`

- [ ] **Step 1: Add failing tests**

Append to test file:

```python
# ── compute_cell_communication() tests ────────────────────────────────────────


def _ha_after_gene_pairs(adata, monkeypatch):
    ha = _setup_ha_no_network(adata, monkeypatch, cell_type_key="cell_type")
    ha._completed_steps.add(STEP_GENE_PAIRS)  # skip actual gene pair computation
    return ha


def test_ccc_standard_marks_step(adata_spatial, monkeypatch):
    ha = _ha_after_gene_pairs(adata_spatial, monkeypatch)
    with patch(
        "scvi.external.harreman._analysis._compute_cell_communication"
    ) as mock_ccc:
        mock_ccc.return_value = None
        ha.compute_cell_communication(mode="standard")
    assert STEP_CCC in ha._completed_steps


def test_ccc_cell_type_mode_calls_ct_function(adata_spatial, monkeypatch):
    ha = _ha_after_gene_pairs(adata_spatial, monkeypatch)
    with patch(
        "scvi.external.harreman._analysis._compute_ct_cell_communication"
    ) as mock_ct:
        mock_ct.return_value = None
        ha.compute_cell_communication(mode="cell_type")
    mock_ct.assert_called_once()


def test_ccc_standard_mode_calls_standard_function(adata_spatial, monkeypatch):
    ha = _ha_after_gene_pairs(adata_spatial, monkeypatch)
    with patch(
        "scvi.external.harreman._analysis._compute_cell_communication"
    ) as mock_std:
        mock_std.return_value = None
        ha.compute_cell_communication(mode="standard")
    mock_std.assert_called_once()


def test_ccc_invalid_mode_raises(adata_spatial, monkeypatch):
    ha = _ha_after_gene_pairs(adata_spatial, monkeypatch)
    with pytest.raises(ValueError, match="mode"):
        ha.compute_cell_communication(mode="invalid")
```

- [ ] **Step 2: Run tests — expect FAIL**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -k "ccc" -v 2>&1 | tail -10
```

- [ ] **Step 3: Add `compute_cell_communication()` to `_analysis.py`**

Add to imports:

```python
from scvi.external.harreman.tools.cell_communication import (
    compute_cell_communication as _compute_cell_communication,
    compute_ct_cell_communication as _compute_ct_cell_communication,
)
```

Add method:

```python
def compute_cell_communication(
    self,
    mode: Literal["standard", "cell_type"] = "standard",
    model: str | None = None,
    n_permutations: int = 1000,
    seed: int = 42,
    test: Literal["parametric", "non-parametric", "both"] = "both",
    mean: Literal["algebraic", "geometric"] = "algebraic",
    device: str | None = None,
    verbose: bool = False,
) -> None:
    """Compute ligand-receptor / metabolite CCC scores.

    Requires ``compute_gene_pairs()`` to have been called.

    Parameters
    ----------
    mode
        ``"standard"`` for cell-type-agnostic analysis;
        ``"cell_type"`` for cell-type-stratified analysis (requires
        ``cell_type_key`` to have been set in ``setup()``).
    model
        Expression normalization model (``"danb"``, ``"normal"``, etc.).
    n_permutations
        Number of permutations for non-parametric test.
    seed
        Random seed for reproducibility.
    test
        Statistical test to run.
    mean
        Averaging method for multi-gene interactions.
    device
        PyTorch device string. ``None`` uses auto-selection.
    verbose
        Print progress.
    """
    self._require(STEP_CCC)
    if mode not in ("standard", "cell_type"):
        raise ValueError(f"mode must be 'standard' or 'cell_type', got '{mode}'")

    kwargs = dict(
        adata=self._adata,
        layer_key_p_test=self._layer_key,
        layer_key_np_test=self._layer_key,
        model=model,
        M=n_permutations,
        seed=seed,
        test=test,
        mean=mean,
        device=device,
        verbose=verbose,
    )
    if mode == "cell_type":
        ct_key = self._adata.uns.get("cell_type_key")
        _compute_ct_cell_communication(**kwargs, cell_type_key=ct_key)
    else:
        _compute_cell_communication(**kwargs)

    self._completed_steps.add(STEP_CCC)
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -k "ccc" -v 2>&1 | tail -15
```

Expected: 4 PASSED

---

## Task 8: `compute_interacting_cell_scores()` + `select_significant_interactions()` + tests

**Files:**
- Modify: `src/scvi/external/harreman/_analysis.py`
- Modify: `tests/external/harreman/test_harreman_analysis.py`

- [ ] **Step 1: Add failing tests**

Append to test file:

```python
# ── compute_interacting_cell_scores() + select_significant_interactions() ──────


def _ha_after_ccc(adata, monkeypatch):
    ha = _ha_after_gene_pairs(adata, monkeypatch)
    ha._completed_steps.add(STEP_CCC)
    return ha


def test_ics_standard_marks_step(adata_spatial, monkeypatch):
    ha = _ha_after_ccc(adata_spatial, monkeypatch)
    with patch(
        "scvi.external.harreman._analysis._compute_interacting_cell_scores"
    ) as mock_ics:
        mock_ics.return_value = None
        ha.compute_interacting_cell_scores(mode="standard")
    assert STEP_ICS in ha._completed_steps


def test_ics_cell_type_mode_calls_ct_function(adata_spatial, monkeypatch):
    ha = _ha_after_ccc(adata_spatial, monkeypatch)
    with patch(
        "scvi.external.harreman._analysis._compute_ct_interacting_cell_scores"
    ) as mock_ct:
        mock_ct.return_value = None
        ha.compute_interacting_cell_scores(mode="cell_type")
    mock_ct.assert_called_once()


def test_select_significant_marks_step(adata_spatial, monkeypatch):
    ha = _ha_after_ccc(adata_spatial, monkeypatch)
    with patch(
        "scvi.external.harreman._analysis._select_significant_interactions"
    ) as mock_sig:
        mock_sig.return_value = None
        ha.select_significant_interactions()
    assert STEP_SIG in ha._completed_steps


def test_select_significant_passes_threshold(adata_spatial, monkeypatch):
    ha = _ha_after_ccc(adata_spatial, monkeypatch)
    with patch(
        "scvi.external.harreman._analysis._select_significant_interactions"
    ) as mock_sig:
        mock_sig.return_value = None
        ha.select_significant_interactions(fdr_threshold=0.01)
    _, kwargs = mock_sig.call_args
    assert kwargs.get("threshold") == 0.01
```

- [ ] **Step 2: Run tests — expect FAIL**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -k "ics or significant" -v 2>&1 | tail -10
```

- [ ] **Step 3: Add both methods to `_analysis.py`**

Add to imports:

```python
from scvi.external.harreman.tools.cell_communication import (
    compute_interacting_cell_scores as _compute_interacting_cell_scores,
    compute_ct_interacting_cell_scores as _compute_ct_interacting_cell_scores,
    select_significant_interactions as _select_significant_interactions,
)
```

Add methods:

```python
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
        Use FDR (``True``) or raw p-values (``False``).
    ct_aware
        Whether to filter cell-type-aware results. Defaults to ``True`` if
        ``mode="cell_type"`` was used, ``False`` otherwise.
    """
    self._require(STEP_SIG)
    if ct_aware is None:
        ct_aware = (
            STEP_CCC in self._completed_steps
            and self._adata.uns.get("cell_type_key") is not None
        )
    _select_significant_interactions(
        adata=self._adata,
        ct_aware=ct_aware,
        test=test,
        use_FDR=use_fdr,
        threshold=fdr_threshold,
    )
    self._completed_steps.add(STEP_SIG)
```

- [ ] **Step 4: Run tests — expect PASS**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -k "ics or significant" -v 2>&1 | tail -15
```

Expected: 4 PASSED

---

## Task 9: Update `__init__.py` + run all unit tests

**Files:**
- Modify: `src/scvi/external/harreman/__init__.py`

- [ ] **Step 1: Update `__init__.py`**

Read current content first, then update:

```python
# src/scvi/external/harreman/__init__.py
from . import datasets as ds
from . import hotspot as hs
from . import preprocessing as pp
from . import tools as tl
from ._analysis import HarremanAnalysis

__all__ = ["ds", "hs", "pp", "tl", "HarremanAnalysis"]
```

- [ ] **Step 2: Verify top-level import**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -c "from scvi.external.harreman import HarremanAnalysis; print(HarremanAnalysis)"
```

Expected: `<class 'scvi.external.harreman._analysis.HarremanAnalysis'>`

- [ ] **Step 3: Run all unit tests**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -v -m "not optional" 2>&1 | tail -30
```

Expected: All unit tests PASS. Count should be ≥ 20 tests.

- [ ] **Step 4: Verify existing harreman tests still pass**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman.py -v 2>&1 | tail -20
```

Expected: All existing tests still pass (no regressions).

---

## Task 10: Integration tests — RESOLVI, SCVIVA, DestVI

**Files:**
- Modify: `tests/external/harreman/test_harreman_analysis.py` (add integration section)

These tests train real models (2 epochs on synthetic data) and verify `HarremanAnalysis.__init__` correctly extracts model outputs. Marked `@pytest.mark.optional` so they don't block CI.

- [ ] **Step 1: Add integration test fixtures and tests**

Append to test file:

```python
# ── Integration tests (require model training) ────────────────────────────────
# Run with: pytest -m optional tests/external/harreman/test_harreman_analysis.py


@pytest.fixture(scope="session")
def destvi_model_and_adata():
    """Train minimal DestVI model on synthetic data."""
    import numpy as np
    from scvi.data import synthetic_iid
    from scvi.model import CondSCVI, DestVI

    n_labels = 3
    dataset = synthetic_iid(n_labels=n_labels, n_genes=50)
    dataset.obs["overclustering_vamp"] = list(range(dataset.n_obs))
    CondSCVI.setup_anndata(dataset, labels_key="labels")
    sc_model = CondSCVI(dataset, n_latent=2, n_layers=1, prior="normal")
    sc_model.train(1, train_size=1)
    del dataset.obs["overclustering_vamp"]
    DestVI.setup_anndata(dataset, layer=None)
    spatial_model = DestVI.from_rna_model(dataset, sc_model, amortization="latent")
    spatial_model.train(max_epochs=1)
    dataset.obsm["spatial"] = np.random.default_rng(0).random((dataset.n_obs, 2)) * 100
    return spatial_model, dataset


@pytest.fixture(scope="session")
def resolvi_model_and_adata():
    """Train minimal RESOLVI model on synthetic spatial data."""
    import numpy as np
    from scvi.data import synthetic_iid
    from scvi.external import RESOLVI

    adata = synthetic_iid(generate_coordinates=True, n_regions=5, n_genes=50)
    adata.obsm["X_spatial"] = adata.obsm["coordinates"]
    RESOLVI.setup_anndata(adata)
    model = RESOLVI(adata)
    model.train(max_epochs=2)
    return model, adata


@pytest.fixture(scope="session")
def scviva_model_and_adata():
    """Train minimal SCVIVA model on synthetic spatial data."""
    import numpy as np
    from scvi.data import synthetic_iid
    from scvi.external import SCVIVA

    adata = synthetic_iid(
        batch_size=64,
        n_genes=50,
        n_proteins=0,
        n_regions=0,
        n_batches=2,
        n_labels=3,
        generate_coordinates=True,
        sparse_format=None,
    )
    adata.layers["counts"] = adata.X.copy()
    adata.obsm["qz1_m"] = np.random.default_rng(1).normal(size=(adata.n_obs, 10))
    SCVIVA.preprocessing_anndata(
        adata,
        k_nn=5,
        sample_key="batch",
        labels_key="labels",
        cell_coordinates_key="coordinates",
        expression_embedding_key="qz1_m",
        expression_embedding_niche_key="qz1_m_niche_ct",
        niche_composition_key="neighborhood_composition",
        niche_indexes_key="niche_indexes",
        niche_distances_key="niche_distances",
    )
    SCVIVA.setup_anndata(
        adata,
        layer="counts",
        batch_key="batch",
        sample_key="batch",
        labels_key="labels",
        cell_coordinates_key="coordinates",
        expression_embedding_key="qz1_m",
        expression_embedding_niche_key="qz1_m_niche_ct",
        niche_composition_key="neighborhood_composition",
        niche_indexes_key="niche_indexes",
        niche_distances_key="niche_distances",
    )
    model = SCVIVA(adata, prior_mixture=False, semisupervised=False)
    model.train(max_epochs=2, accelerator="cpu")
    return model, adata


@pytest.mark.optional
def test_integration_destvi_attaches_ct_layers(destvi_model_and_adata):
    """DestVI: verify each cell-type layer is attached to adata."""
    model, adata = destvi_model_and_adata
    ha = HarremanAnalysis(adata.copy(), model=model)
    proportions = model.get_proportions()
    for ct in proportions.columns:
        assert ct in ha.adata.layers, f"Missing layer for cell type '{ct}'"


@pytest.mark.optional
def test_integration_destvi_is_deconvolved(destvi_model_and_adata):
    model, adata = destvi_model_and_adata
    ha = HarremanAnalysis(adata.copy(), model=model)
    assert ha.is_deconvolved is True


@pytest.mark.optional
def test_integration_destvi_layer_shape(destvi_model_and_adata):
    """Cell-type layers must have same shape as adata.X."""
    model, adata = destvi_model_and_adata
    adata_copy = adata.copy()
    ha = HarremanAnalysis(adata_copy, model=model)
    proportions = model.get_proportions()
    for ct in proportions.columns:
        assert ha.adata.layers[ct].shape == adata.shape, f"Layer '{ct}' shape mismatch"


@pytest.mark.optional
def test_integration_resolvi_attaches_denoised_layer(resolvi_model_and_adata):
    """RESOLVI: verify denoised layer attached."""
    from scvi.external.harreman._constants import HARREMAN_DENOISED_LAYER

    model, adata = resolvi_model_and_adata
    ha = HarremanAnalysis(adata.copy(), model=model)
    assert HARREMAN_DENOISED_LAYER in ha.adata.layers


@pytest.mark.optional
def test_integration_resolvi_not_deconvolved(resolvi_model_and_adata):
    model, adata = resolvi_model_and_adata
    ha = HarremanAnalysis(adata.copy(), model=model)
    assert ha.is_deconvolved is False


@pytest.mark.optional
def test_integration_resolvi_denoised_shape(resolvi_model_and_adata):
    from scvi.external.harreman._constants import HARREMAN_DENOISED_LAYER

    model, adata = resolvi_model_and_adata
    adata_copy = adata.copy()
    ha = HarremanAnalysis(adata_copy, model=model)
    assert ha.adata.layers[HARREMAN_DENOISED_LAYER].shape == adata.shape


@pytest.mark.optional
def test_integration_scviva_attaches_latent(scviva_model_and_adata):
    """SCVIVA: verify latent representation attached to obsm."""
    from scvi.external.harreman._constants import HARREMAN_LATENT_OBSM

    model, adata = scviva_model_and_adata
    ha = HarremanAnalysis(adata.copy(), model=model)
    assert HARREMAN_LATENT_OBSM in ha.adata.obsm


@pytest.mark.optional
def test_integration_scviva_not_deconvolved(scviva_model_and_adata):
    model, adata = scviva_model_and_adata
    ha = HarremanAnalysis(adata.copy(), model=model)
    assert ha.is_deconvolved is False


@pytest.mark.optional
def test_integration_scviva_latent_shape(scviva_model_and_adata):
    from scvi.external.harreman._constants import HARREMAN_LATENT_OBSM

    model, adata = scviva_model_and_adata
    adata_copy = adata.copy()
    ha = HarremanAnalysis(adata_copy, model=model)
    latent = ha.adata.obsm[HARREMAN_LATENT_OBSM]
    assert latent.shape[0] == adata.n_obs
    assert latent.ndim == 2
```

- [ ] **Step 2: Run integration tests**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -m optional -v 2>&1 | tail -30
```

Expected: 9 PASSED (DestVI ×3, RESOLVI ×3, SCVIVA ×3)

- [ ] **Step 3: Confirm unit tests still pass**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/test_harreman_analysis.py -v 2>&1 | tail -20
```

Expected: All PASSED. Print total count.

---

## Task 11: Ruff / linting pass

**Files:**
- Modify: `src/scvi/external/harreman/_analysis.py`, `_results.py`, `_constants.py`

- [ ] **Step 1: Run ruff check on new files**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m ruff check src/scvi/external/harreman/_analysis.py src/scvi/external/harreman/_results.py src/scvi/external/harreman/_constants.py 2>&1
```

Fix any reported violations before moving on.

- [ ] **Step 2: Run ruff format**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m ruff format src/scvi/external/harreman/_analysis.py src/scvi/external/harreman/_results.py src/scvi/external/harreman/_constants.py tests/external/harreman/test_harreman_analysis.py
```

- [ ] **Step 3: Re-run all tests to confirm no regressions**

```bash
cd /Users/orikr/PycharmProjects/harreman
python -m pytest tests/external/harreman/ -v 2>&1 | tail -30
```

Expected: All unit tests PASS, existing tests PASS.

---

## Self-Review Checklist

- [x] **`_constants.py`** — Task 1 covers all constants used in later tasks
- [x] **`_results.py`** — Task 2 covers `HarremanResults.from_uns`, no other constructors used
- [x] **Init/properties** — Task 3 covers TypeError, ValueError (unsupported model), RuntimeError (step ordering), all properties
- [x] **`setup()`** — Task 4 covers KNN graph build + database extraction (mocked). Params written to `adata.uns`
- [x] **`filter_genes()`** — Task 5 covers prerequisite check + delegation
- [x] **`compute_gene_pairs()`** — Task 6 covers prerequisite check + layer_key passthrough
- [x] **`compute_cell_communication()`** — Task 7 covers both modes + invalid mode
- [x] **`compute_interacting_cell_scores()`** — Task 8 covers both modes
- [x] **`select_significant_interactions()`** — Task 8 covers threshold passthrough
- [x] **`__init__.py` export** — Task 9
- [x] **DestVI integration** — Task 10: layers attached, is_deconvolved=True, shapes correct
- [x] **RESOLVI integration** — Task 10: denoised layer attached, is_deconvolved=False, shape correct
- [x] **SCVIVA integration** — Task 10: latent obsm attached, is_deconvolved=False, shape correct
- [x] **Linting** — Task 11

**Type consistency check:**
- `HARREMAN_DENOISED_LAYER` defined in `_constants.py` Task 1, used in `_analysis.py` Task 3, referenced in integration tests Task 10 ✓
- `HARREMAN_LATENT_OBSM` same path ✓
- `_extract_interaction_db` imported at module level in `_analysis.py`, monkeypatched by name in tests ✓
- `HarremanResults.from_uns` defined in Task 2, called in `results` property in Task 3 ✓

**Spec gaps:** Spec mentioned `database_path` parameter — **corrected** in plan: HarremanDB is downloaded via pooch automatically. `setup()` uses `species` + `database` (type) instead.
