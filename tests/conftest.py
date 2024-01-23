import shutil
from distutils.dir_util import copy_tree

import pytest

import scvi
from tests.data.utils import generic_setup_adata_manager


def pytest_addoption(parser):
    """Docstring for pytest_addoption."""
    parser.addoption(
        "--slow",
        action="store_true",
        default=False,
        help="Run slow tests.",
    )
    parser.addoption(
        "--private",
        action="store_true",
        default=False,
        help="Run tests that are private (sensitive credentials or code).",
    )
    parser.addoption(
        "--unit",
        action="store_true",
        default=False,
        help="Run unit tests.",
    )
    parser.addoption(
        "--integration",
        action="store_true",
        default=False,
        help="Run integration tests.",
    )
    parser.addoption(
        "--distributed",
        action="store_true",
        default=False,
        help="Run distributed tests.",
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
    config.addinivalue_line("markers", "slow: mark test as slow.")
    config.addinivalue_line("markers", "private: mark test as private.")
    config.addinivalue_line("markers", "unit: mark test as unit.")
    config.addinivalue_line("markers", "integration: mark test as integration.")


def pytest_collection_modifyitems(config, items):
    """Docstring for pytest_collection_modifyitems."""
    run_slow = config.getoption("--slow")
    skip_slow = pytest.mark.skip(reason="need --slow option to run")
    for item in items:
        # All tests marked with `pytest.mark.slow` get skipped unless `--slow` passed
        if not run_slow and ("slow" in item.keywords):
            item.add_marker(skip_slow)

    run_private = config.getoption("--private")
    skip_private = pytest.mark.skip(reason="need --private option to run")
    for item in items:
        # All tests marked with `pytest.mark.private` get skipped unless
        # `--private` passed
        if not run_private and ("private" in item.keywords):
            item.add_marker(skip_private)
        # Skip all tests not marked with `pytest.mark.private` if `--private` passed
        elif run_private and ("private" not in item.keywords):
            item.add_marker(skip_private)

    run_unit = config.getoption("--unit")
    skip_unit = pytest.mark.skip(reason="need --unit option to run")
    for item in items:
        # All tests marked with `pytest.mark.unit` get skipped unless `--unit` passed
        if not run_unit and ("unit" in item.keywords):
            item.add_marker(skip_unit)

    run_integration = config.getoption("--integration")
    skip_integration = pytest.mark.skip(reason="need --integration option to run")
    for item in items:
        # All tests marked with `pytest.mark.integration` get skipped unless
        # `--integration` passed
        if not run_integration and ("integration" in item.keywords):
            item.add_marker(skip_integration)

    run_distributed = config.getoption("--distributed")
    skip_distributed = pytest.mark.skip(reason="need --distributed option to run")
    for item in items:
        # All tests marked with `pytest.mark.distributed` get skipped unless
        # `--distributed` passed
        if not run_distributed and ("distributed" in item.keywords):
            item.add_marker(skip_distributed)
        # Skip all tests not marked with `pytest.mark.distributed` if `--distributed` passed
        elif run_distributed and ("distributed" not in item.keywords):
            item.add_marker(skip_distributed)


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


@pytest.fixture(scope="session")
def mock_contrastive_adata():
    """Synthetic contrastive adata."""
    adata = scvi.data.synthetic_iid(n_batches=2)
    adata = adata[:-3, :]  # Unequal technical batch sizes.
    adata.layers["raw_counts"] = adata.X.copy()
    return adata


@pytest.fixture(scope="session")
def mock_contrastive_adata_manager(mock_contrastive_adata):
    """Anndata manager for synthetic contrastive data."""
    return generic_setup_adata_manager(
        adata=mock_contrastive_adata,
        batch_key="batch",
        labels_key="labels",
        layer="raw_counts",
    )


@pytest.fixture(scope="session")
def mock_background_indices(mock_contrastive_adata):
    """Indices for background data in ``mock_contrastive_adata_manager``."""
    adata = mock_contrastive_adata
    return adata.obs.index[(adata.obs["batch"] == "batch_0")].astype(int).tolist()


@pytest.fixture(scope="session")
def mock_target_indices(mock_contrastive_adata):
    """Indices for target data in ``mock_contrastive_adata_manager``."""
    adata = mock_contrastive_adata
    return adata.obs.index[(adata.obs["batch"] == "batch_1")].astype(int).tolist()
