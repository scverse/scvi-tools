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
    parser.addoption(
        "--optional",
        action="store_true",
        default=False,
        help="Run tests that are optional.",
    )


def pytest_configure(config):
    """Docstring for pytest_configure."""
    config.addinivalue_line("markers", "optional: mark test as optional.")


def pytest_collection_modifyitems(config, items):
    """Docstring for pytest_collection_modifyitems."""
    run_internet = config.getoption("--internet-tests")
    skip_internet = pytest.mark.skip(reason="need --internet-tests option to run")
    for item in items:
        # All tests marked with `pytest.mark.internet` get skipped unless
        # `--internet-tests` passed
        if not run_internet and ("internet" in item.keywords):
            item.add_marker(skip_internet)

    run_optional = config.getoption("--optional")
    skip_optional = pytest.mark.skip(reason="need --optional option to run")
    for item in items:
        # All tests marked with `pytest.mark.optional` get skipped unless
        # `--optional` passed
        if not run_optional and ("optional" in item.keywords):
            item.add_marker(skip_optional)


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
