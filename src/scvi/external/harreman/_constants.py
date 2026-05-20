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
