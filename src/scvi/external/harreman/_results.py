from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

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
    from typing import Any

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
