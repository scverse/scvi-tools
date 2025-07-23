from importlib.metadata import metadata

import pytest

import scvi
from scvi.hub._url import validate_colab_notebook, validate_url

info = metadata("scvi-tools")


@pytest.mark.internet
def test_validate_homepage_url():
    valid = "https://scvi-tools.org/"
    assert validate_url(valid)

    invalid = "scvi-tools.org"
    assert not validate_url(invalid)


@pytest.mark.internet
def test_validate_hf_url():
    valid = "https://huggingface.co/scvi-tools/"
    assert validate_url(valid)

    invalid = "huggingface.co/scvi-tools"
    assert not validate_url(invalid)


def bump_patch_version(version: str) -> str:
    parts = version.split(".")
    parts[-1] = str(int(parts[-1]) + 1)
    return ".".join(parts)


@pytest.mark.internet
def test_validate_colab_url():
    colab_url1 = (
        "https://colab.research.google.com/github/scverse/scvi-tutorials/blob/"
        + scvi.__version__
        + "/quick_start/api_overview.ipynb"
    )
    assert validate_colab_notebook(colab_url1)

    colab_url2 = (
        "https://colab.research.google.com/github/scverse/scvi-tutorials/blob/"
        + bump_patch_version(scvi.__version__)
        + "/quick_start/api_overview.ipynb"
    )
    assert not validate_colab_notebook(colab_url2)
