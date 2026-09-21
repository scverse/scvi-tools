import ast
import pathlib
import warnings

import numpy as np
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
        tree = ast.parse(path.read_text(encoding="utf-8"))
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


def _dataloader_arg_warnings(record):
    """Messages complaining that an argument is redundant next to a custom dataloader.

    Matched loosely on purpose so that the assertions below bind to the behaviour
    rather than to one exact phrasing.
    """
    messages = []
    for w in record:
        message = str(w.message)
        lowered = message.lower()
        if "dataloader" in lowered and ("redundant" in lowered or "ignored" in lowered):
            messages.append(message)
    return messages


@pytest.fixture(scope="module")
def trained_scvi_model():
    adata = scvi.data.synthetic_iid()
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    model = scvi.model.SCVI(adata, n_latent=5)
    model.train(max_epochs=1, train_size=1.0, accelerator="cpu")
    return model


@pytest.mark.parametrize(
    "method",
    [
        "get_normalized_expression",
        "get_likelihood_parameters",
        "get_latent_representation",
        "get_elbo",
    ],
)
def test_dataloader_does_not_warn_for_untouched_args(trained_scvi_model, method):
    """Passing only a ``dataloader`` must not warn about arguments the caller never set.

    ``n_samples`` defaults to a non-``None`` value, so a ``value is not None`` check
    reported it as redundant on every call.
    """
    model = trained_scvi_model
    dataloader = model._make_data_loader(model.adata)

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        getattr(model, method)(dataloader=dataloader)

    assert _dataloader_arg_warnings(record) == []


def test_dataloader_warns_by_argument_name(trained_scvi_model):
    """An explicitly set argument is reported by name, not by value."""
    model = trained_scvi_model
    dataloader = model._make_data_loader(model.adata)

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        model.get_normalized_expression(dataloader=dataloader, batch_size=128)

    messages = _dataloader_arg_warnings(record)
    assert len(messages) == 1
    assert "`batch_size`" in messages[0]
    assert "128" not in messages[0]


def test_dataloader_warns_for_array_indices(trained_scvi_model):
    """``indices`` may be an array; the check must not raise on ambiguous truth values."""
    model = trained_scvi_model
    dataloader = model._make_data_loader(model.adata)

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        model.get_latent_representation(dataloader=dataloader, indices=np.arange(5))

    messages = _dataloader_arg_warnings(record)
    assert len(messages) == 1
    assert "`indices`" in messages[0]
