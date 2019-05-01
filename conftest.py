import pytest
from distutils.dir_util import copy_tree
import shutil


def pytest_addoption(parser):
    parser.addoption("--model_fit", action="store_true",
                     help="run the tests only in case of that command line (marked with marker @no_cmd)")


@pytest.fixture(scope="session")
def save_path(tmpdir_factory):
    dir = tmpdir_factory.mktemp("temp_data", numbered=False)
    path = str(dir)
    copy_tree('tests/data', path)
    yield path + '/'
    shutil.rmtree(str(tmpdir_factory.getbasetemp()))


def pytest_runtest_setup(item):
    if 'model_fit' in item.keywords and not item.config.getoption("--model_fit"):
        pytest.skip("need --model_fit option to run this test")
