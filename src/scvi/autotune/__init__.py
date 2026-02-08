from __future__ import annotations

import os
from typing import TYPE_CHECKING

from scvi.utils import error_on_missing_dependencies

if TYPE_CHECKING:
    from ._experiment import AutotuneExperiment, ScibTuneReportCheckpointCallback
    from ._tune import run_autotune

__all__ = ["AutotuneExperiment", "ScibTuneReportCheckpointCallback", "run_autotune"]

_DOCS_BUILD = os.environ.get("READTHEDOCS") == "True"

if _DOCS_BUILD:
    _AUTOTUNE_IMPORT_ERROR = ModuleNotFoundError(
        "Autotune requires optional dependencies; install scvi-tools[autotune]."
    )

    class _MissingAutotuneDependency:
        """Autotune is unavailable because dependencies failed to import."""

        def __init__(self, *args, **kwargs):
            raise _AUTOTUNE_IMPORT_ERROR

    class AutotuneExperiment(_MissingAutotuneDependency):
        """Autotune is unavailable; see class docstring."""

    class ScibTuneReportCheckpointCallback(_MissingAutotuneDependency):
        """Autotune is unavailable; see class docstring."""

    def run_autotune(*args, **kwargs):
        """Autotune is unavailable; install scvi-tools[autotune]."""
        raise _AUTOTUNE_IMPORT_ERROR
else:
    try:
        error_on_missing_dependencies("hyperopt", "ray.tune")
        from ._experiment import AutotuneExperiment, ScibTuneReportCheckpointCallback
        from ._tune import run_autotune
    except ImportError as _autotune_import_error:  # pragma: no cover - docs/optional deps
        _AUTOTUNE_IMPORT_ERROR = _autotune_import_error

        class _MissingAutotuneDependency:
            """Autotune is unavailable because dependencies failed to import."""

            def __init__(self, *args, **kwargs):
                raise _AUTOTUNE_IMPORT_ERROR

        class AutotuneExperiment(_MissingAutotuneDependency):
            """Autotune is unavailable; see class docstring."""

        class ScibTuneReportCheckpointCallback(_MissingAutotuneDependency):
            """Autotune is unavailable; see class docstring."""

        def run_autotune(*args, **kwargs):
            """Autotune is unavailable; install scvi-tools[autotune]."""
            raise _AUTOTUNE_IMPORT_ERROR
