from __future__ import annotations

import pytest


# TODO: use fixtures for models, anndata that will be used in tests
@pytest.fixture(scope="session")
def adata():
    pass


def model():
    pass


@pytest.mark.parametrize
def test_get_aggregated_posterior():
    pass


# TODO: use parametrize to test multiple values for the parameters
@pytest.mark.parametrize

# TODO: define the actual test
def test_differential_abundance():
    pass


# TODO: repeat the above 2 steps as many times as needed (may actually only need one test)
