from scvi.utils import error_on_missing_dependencies

from ._decorators import auto_move_data

error_on_missing_dependencies("pyro")

from pyro.infer import Predictive  # noqa: E402


class AutoMoveDataPredictive(Predictive):
    """Auto move data."""

    @auto_move_data
    def forward(self, *args, **kwargs):
        return super().forward(*args, **kwargs)
