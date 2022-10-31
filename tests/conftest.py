import shutil
from distutils.dir_util import copy_tree

import pytest

import scvi


def pytest_addoption(parser):
    """Docstring for pytest_addoption."""
    parser.addoption(
        "--model_fit",
        action="store_true",
        default=False,
        dest="model_fit",
        help="Option to run full training model for test_model_fit",
    )
    parser.addoption(
        "--internet-tests",
        action="store_true",
        default=False,
        help="Run tests that retrieve stuff from the internet. This increases test time.",
    )


def pytest_collection_modifyitems(config, items):
    """Docstring for pytest_collection_modifyitems."""
    run_internet = config.getoption("--internet-tests")
    skip_internet = pytest.mark.skip(reason="need --internet-tests option to run")
    for item in items:
        # All tests marked with `pytest.mark.internet` get skipped unless
        # `--internet-tests` passed
        if not run_internet and ("internet" in item.keywords):
            item.add_marker(skip_internet)


@pytest.fixture(scope="session")
def save_path(tmpdir_factory):
    """Docstring for save_path."""
    dir = tmpdir_factory.mktemp("temp_data", numbered=False)
    path = str(dir)
    copy_tree("tests/data", path)
    yield path + "/"
    shutil.rmtree(str(tmpdir_factory.getbasetemp()))


@pytest.fixture(scope="session")
def synthetic_adata():
    """Docstring for model_fit."""
    return scvi.data.synthetic_iid()


@pytest.fixture(scope="session")
def model_fit(request):
    """Docstring for model_fit."""
    return request.config.getoption("--model_fit")
