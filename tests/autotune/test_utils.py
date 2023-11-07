import pytest

import scvi


@pytest.mark.optional
def test_in_notebook():
    assert not scvi.autotune._utils.in_notebook()
