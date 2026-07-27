import inspect
import re

import pytest

import scvi.model
from scvi.external import MRVI, RESOLVI, ContrastiveVI

# Public model classes whose method docstrings should be fully rendered, i.e.
# contain no leftover ``%(...)s`` docstring-template placeholders. The external
# models listed here define methods that share the ``de_dsp`` templated
# parameter descriptions.
_MODEL_CLASSES = [
    getattr(scvi.model, name)
    for name in scvi.model.__all__
    if inspect.isclass(getattr(scvi.model, name))
] + [ContrastiveVI, MRVI, RESOLVI]

_PLACEHOLDER = re.compile(r"%\([A-Za-z0-9_]+\)s")


@pytest.mark.parametrize("model_cls", _MODEL_CLASSES, ids=lambda c: c.__name__)
def test_no_unsubstituted_docstring_placeholders(model_cls):
    """Public methods must not expose raw ``%(...)s`` docstring placeholders.

    A method whose docstring references a ``de_dsp`` placeholder (e.g.
    ``%(de_silent)s``) but is missing the ``@de_dsp.dedent`` decorator renders
    the literal placeholder to users. This guards against that regression.
    """
    offenders = []
    for name, member in inspect.getmembers(model_cls, predicate=inspect.isfunction):
        if name.startswith("_"):
            continue
        doc = member.__doc__ or ""
        leftover = _PLACEHOLDER.findall(doc)
        if leftover:
            offenders.append(f"{model_cls.__name__}.{name}: {sorted(set(leftover))}")
    assert not offenders, "Unsubstituted docstring placeholders found:\n" + "\n".join(offenders)
