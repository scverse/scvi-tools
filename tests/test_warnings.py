import ast
import pathlib

import pytest

import scvi

_SRC = pathlib.Path(scvi.__file__).parent


def _find_noop_warning_statements():
    """Return ``file:lineno`` for every bare ``*Warning(...)`` statement in the source.

    A statement such as ``Warning("msg")`` merely constructs a warning object and
    immediately discards it, so the intended user-facing warning never fires. The
    correct form is ``warnings.warn("msg", ...)``.
    """
    offenders = []
    for path in _SRC.rglob("*.py"):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
                func = node.value.func
                name = getattr(func, "id", None) or getattr(func, "attr", None)
                if name is not None and name.endswith("Warning"):
                    offenders.append(f"{path.relative_to(_SRC)}:{node.lineno}")
    return offenders


def test_no_bare_warning_statements():
    offenders = _find_noop_warning_statements()
    assert not offenders, (
        "Found `Warning(...)` used as a no-op statement (should be `warnings.warn(...)`):\n"
        + "\n".join(offenders)
    )


def test_get_ranked_features_warns_on_missing_attrs():
    adata = scvi.data.synthetic_iid()
    scvi.model.SCANVI.setup_anndata(
        adata, labels_key="labels", unlabeled_category="label_0", batch_key="batch"
    )
    model = scvi.model.SCANVI(adata)
    with pytest.warns(UserWarning, match="Missing Attributions matrix"):
        result = model.get_ranked_features()
    assert result is None
