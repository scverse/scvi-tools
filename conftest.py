import pytest
from distutils.dir_util import copy_tree
import shutil


def pytest_addoption(parser):
    parser.addoption(
        "--model_fit",
        action="store_true",
        default=False,
        dest="model_fit",
        help="Option to run full training model for test_model_fit",
    )


@pytest.fixture(scope="session")
def save_path(tmpdir_factory):
    dir = tmpdir_factory.mktemp("temp_data", numbered=False)
    path = str(dir)
    copy_tree("tests/data", path)
    yield path + "/"
    shutil.rmtree(str(tmpdir_factory.getbasetemp()))


@pytest.fixture(scope="session")
def model_fit(request):
    return request.config.getoption("--model_fit")
