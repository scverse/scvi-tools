from scvi.utils import in_notebook


def test_in_notebook():
    assert not in_notebook()
