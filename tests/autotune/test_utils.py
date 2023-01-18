import scvi


def test_in_notebook():
    assert not scvi.autotune.in_notebook()
