---
name: harreman-block-diagram-design
description: Design spec for the Harreman HTML block diagram showing HarremanAnalysis at center connected to spatial models and submodules
metadata:
  type: project
---

# Harreman Block Diagram Design

**Date:** 2026-06-01
**Status:** Approved
**Output:** `docs/architecture/harreman-block-diagram.html`

---

## Layout: Three-tier grid

### Tier 1 — Spatial model inputs (top row, 4 cards)

| Card | Color token | Key detail |
|------|-------------|------------|
| DestVI | `--external` (amber) | `get_proportions()` → cell-type count layers; enables `mode="cell_type"` |
| ResolVI | `--external` (amber) | `get_normalized_expression()` → `layers["denoised"]` |
| SCVIVA | `--external` (amber) | `get_latent_representation()` → latent coords for KNN |
| Standalone | `--upstream` (gray) | pre-prepared AnnData, no model needed |

All four feed into HarremanAnalysis via a downward arrow ("feeds into").

### Tier 2 — HarremanAnalysis (full-width hub card, `--foundation` teal)

Central orchestrator class. Chips:
`setup()` · `filter_genes()` · `compute_gene_pairs()` · `compute_cell_communication()` · `compute_interacting_cell_scores()` · `select_significant_interactions()` · `ha.hs` · `ha.tl` · `ha.pl`

State stored in `adata.uns["harreman"]`.

### Tier 3 — Submodule groups (3-column grid)

| Column | Submodules | Colors |
|--------|------------|--------|
| Left | `preprocessing` + `datasets` | `--upstream` gray + `--data` yellow |
| Center | `hotspot` | `--mixin` green |
| Right | `tools` + `plots` | `--core` blue + `--imaging` purple |

HarremanAnalysis connects down to each column with "wraps" dashed arrows.

### Data flow rail (bottom)

Five steps, left-to-right arrows:
1. **Inputs** — AnnData / spatial coords / HarremanDB
2. **Register** — `preprocessing.setup_anndata()`
3. **Autocorrelation** — `hotspot.compute_local_autocorrelation()` + `compute_local_correlation()`
4. **CCC inference** — `tools.compute_cell_communication()` + `compute_interacting_cell_scores()`
5. **Outputs** — `adata.uns["harreman"]` · plots · `HarremanResults`

---

## Color scheme

Reuses the same CSS custom-property palette as `scviva-tools-block-diagram.html`:

| Token | Hex | Semantic use |
|-------|-----|--------------|
| `--upstream` | `#cfd8df` | upstream/standalone inputs, preprocessing |
| `--foundation` | `#8fd3c7` | HarremanAnalysis hub |
| `--mixin` | `#b8d878` | hotspot (composed statistical layer) |
| `--core` | `#8fb7e3` | tools |
| `--external` | `#efbd73` | spatial model inputs (DestVI, ResolVI, SCVIVA) |
| `--imaging` | `#c4aadd` | plots |
| `--data` | `#f0dd86` | data flow rail + datasets |

## Arrow types

- Solid downward arrow → "feeds into" (models → HarremanAnalysis)
- Dashed downward arrow → "wraps" (HarremanAnalysis → submodules)
- Solid horizontal arrow → data flow (rail steps)
