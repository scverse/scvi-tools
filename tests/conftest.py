import shutil
from distutils.dir_util import copy_tree

import pytest


def pytest_addoption(parser):
    """Docstring for pytest_addoption."""
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
    parser.addoption(
        "--accelerator",
        action="store",
        default="cpu",
        help="Option to specify which accelerator to use for tests.",
    )
    parser.addoption(
        "--devices",
        action="store",
        default="auto",
        help="Option to specify which devices to use for tests.",
    )
    parser.addoption(
        "--seed",
        action="store",
        default=0,
        help="Option to specify which scvi-tools seed to use for tests.",
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
def save_path(tmp_path_factory):
    """Docstring for save_path."""
    dir = tmp_path_factory.mktemp("temp_data", numbered=False)
    path = str(dir)
    copy_tree("tests/test_data", path)
    yield path + "/"
    shutil.rmtree(str(tmp_path_factory.getbasetemp()))


@pytest.fixture(scope="session")
def accelerator(request):
    """Docstring for accelerator."""
    return request.config.getoption("--accelerator")


@pytest.fixture(scope="session")
def devices(request):
    """Docstring for devices."""
    return request.config.getoption("--devices")


@pytest.fixture(autouse=True)
def set_seed(request):
    """Sets the seed for each test."""
    from scvi import settings

    settings.seed = int(request.config.getoption("--seed"))
    yield
