import scvi


def test_in_notebook():
    assert not scvi.autotune._utils.in_notebook()
