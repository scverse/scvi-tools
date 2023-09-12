import shutil
from distutils.dir_util import copy_tree

import pytest

import scvi
from tests.dataset.utils import generic_setup_adata_manager


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


@pytest.fixture(scope="session")
def accelerator(request):
    """Docstring for accelerator."""
    return request.config.getoption("--accelerator")


@pytest.fixture(scope="session")
def devices(request):
    """Docstring for devices."""
    return request.config.getoption("--devices")


@pytest.fixture(scope="session")
def mock_contrastive_adata_manager():
    """adata manager for synthetic contrastive data."""
    adata = scvi.data.synthetic_iid(n_batches=2)
    adata = adata[:-3, :]  # Unequal technical batch sizes.
    adata.layers["raw_counts"] = adata.X.copy()
    return generic_setup_adata_manager(
        adata=adata, batch_key="batch", labels_key="labels", layer="raw_counts"
    )


@pytest.fixture(scope="session")
def mock_background_indices(mock_contrastive_adata_manager):
    """Indices for background data in ``mock_contrastive_adata_manager``."""
    adata = mock_contrastive_adata_manager.adata
    return adata.obs.index[(adata.obs["batch"] == "batch_0")].astype(int).tolist()


@pytest.fixture(scope="session")
def mock_target_indices(mock_contrastive_adata_manager):
    """Indices for target data in ``mock_contrastive_adata_manager``."""
    adata = mock_contrastive_adata_manager.adata
    return adata.obs.index[(adata.obs["batch"] == "batch_1")].astype(int).tolist()
