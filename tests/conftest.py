import shutil

import pytest
from distutils.dir_util import copy_tree

import scvi
from tests.data.utils import generic_setup_adata_manager


def pytest_addoption(parser):
    """Docstring for pytest_addoption."""
    parser.addoption(
        "--internet-tests",
        action="store_true",
        default=False,
        help="Run tests that retrieve stuff from the internet. This increases test time.",
    )
    parser.addoption(
        "--multigpu-tests",
        action="store_true",
        default=False,
        help="Run tests that are designed for multiGPU.",
    )
    parser.addoption(
        "--autotune-tests",
        action="store_true",
        default=False,
        help="Run tests that are designed for Ray Autotune.",
    )
    parser.addoption(
        "--custom-dataloader-tests",
        action="store_true",
        default=False,
        help="Run tests that deals with custom dataloaders. This increases test time.",
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
    parser.addoption(
        "--private",
        action="store_true",
        default=False,
        help="Run tests that are private.",
    )


def pytest_configure(config):
    """Docstring for pytest_configure."""
    config.addinivalue_line("markers", "optional: mark test as optional.")


def pytest_collection_modifyitems(config, items):
    """Docstring for pytest_collection_modifyitems."""
    run_internet = config.getoption("--internet-tests")
    skip_internet = pytest.mark.skip(reason="need --internet-tests option to run")
    skip_non_internet = pytest.mark.skip(reason="test not having a pytest.mark.internet decorator")
    for item in items:
        # All tests marked with `pytest.mark.internet` get skipped unless
        # `--internet-tests` passed
        if not run_internet and ("internet" in item.keywords):
            item.add_marker(skip_internet)
        # Skip all tests not marked with `pytest.mark.internet` if `--internet-tests` passed
        elif run_internet and ("internet" not in item.keywords):
            item.add_marker(skip_non_internet)

    run_custom_dataloader = config.getoption("--custom-dataloader-tests")
    skip_custom_dataloader = pytest.mark.skip(
        reason="need ---custom-dataloader-tests option to run"
    )
    skip_non_custom_dataloader = pytest.mark.skip(
        reason="test not having a pytest.mark.custom_dataloader decorator"
    )
    for item in items:
        # All tests marked with `pytest.mark.custom_dataloader` get skipped unless
        # `--custom_dataloader-tests` passed
        if not run_custom_dataloader and ("dataloader" in item.keywords):
            item.add_marker(skip_custom_dataloader)
        # Skip all tests not marked with `pytest.mark.custom_dataloader`
        # if `--custom-dataloader-tests` passed
        elif run_custom_dataloader and ("dataloader" not in item.keywords):
            item.add_marker(skip_non_custom_dataloader)

    run_optional = config.getoption("--optional")
    skip_optional = pytest.mark.skip(reason="need --optional option to run")
    skip_non_optional = pytest.mark.skip(reason="test not having a pytest.mark.optional decorator")
    for item in items:
        # All tests marked with `pytest.mark.optional` get skipped unless
        # `--optional` passed
        if not run_optional and ("optional" in item.keywords):
            item.add_marker(skip_optional)
        # Skip all tests not marked with `pytest.mark.optional` if `--optional` passed
        elif run_optional and ("optional" not in item.keywords):
            item.add_marker(skip_non_optional)

    run_private = config.getoption("--private")
    skip_private = pytest.mark.skip(reason="need --private option to run")
    skip_non_private = pytest.mark.skip(reason="test not having a pytest.mark.private decorator")
    for item in items:
        # All tests marked with `pytest.mark.private` get skipped unless
        # `--private` passed
        if not run_private and ("private" in item.keywords):
            item.add_marker(skip_private)
        # Skip all tests not marked with `pytest.mark.private` if `--private` passed
        elif run_private and ("private" not in item.keywords):
            item.add_marker(skip_non_private)

    run_multigpu = config.getoption("--multigpu-tests")
    skip_multigpu = pytest.mark.skip(reason="need --multigpu-tests option to run")
    skip_non_multigpu = pytest.mark.skip(reason="test not having a pytest.mark.multigpu decorator")
    for item in items:
        # All tests marked with `pytest.mark.multigpu` get skipped unless
        # `--multigpu-tests` passed
        if not run_multigpu and ("multigpu" in item.keywords):
            item.add_marker(skip_multigpu)
        # Skip all tests not marked with `pytest.mark.multigpu` if `--multigpu-tests` passed
        elif run_multigpu and ("multigpu" not in item.keywords):
            item.add_marker(skip_non_multigpu)

    run_autotune = config.getoption("--autotune-tests")
    skip_autotune = pytest.mark.skip(reason="need --autotune-tests option to run")
    skip_non_autotune = pytest.mark.skip(reason="test not having a pytest.mark.autotune decorator")
    for item in items:
        # All tests marked with `pytest.mark.autotune` get skipped unless
        # `--autotune-tests` passed
        if not run_autotune and ("autotune" in item.keywords):
            item.add_marker(skip_autotune)
        # Skip all tests not marked with `pytest.mark.autotune` if `--autotune-tests` passed
        elif run_autotune and ("autotune" not in item.keywords):
            item.add_marker(skip_non_autotune)


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


@pytest.fixture(autouse=True)
def set_seed(request):
    """Sets the seed for each test."""
    from scvi import settings

    settings.seed = int(request.config.getoption("--seed"))
    return None
