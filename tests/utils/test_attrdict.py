import pytest

from scvi.utils import attrdict


def test_attrdict_new():
    ad = attrdict()
    ad.a = 1
    assert ad.a == 1 and ad["a"] == 1

    ad["b"] = "test"
    assert ad.b == "test" and ad["b"] == "test"

    ad.c = {"d": 2}
    assert isinstance(ad.c, dict) and ad.c["d"] == 2
    with pytest.raises(AttributeError):
        print(ad.c.d)


def test_attrdict_from_dict():
    dct = {"a": 1, "b": 2, "c": {"d": 3}}
    ad = attrdict(dct)
    assert isinstance(ad, attrdict) and not isinstance(ad.c, attrdict)
    assert ad.a == 1 and ad.b == 2 and ad.c["d"] == 3
    assert dict(ad) == dct

    # attrdict creates a copy of dct.
    ad.a = 2
    assert ad.a == 2 and dct["a"] == 1
    ad.c["d"] = 4
    assert ad.c["d"] == 4 and dct["c"]["d"] == 3


def test_attrdict_recursive():
    dct = {"a": 1, "b": 2, "c": {"d": 3}}
    ad = attrdict(dct, recursive=True)
    assert isinstance(ad, attrdict) and isinstance(ad.c, attrdict)
    assert ad.a == 1 and ad.b == 2 and ad.c.d == 3
    assert dict(ad) == dct

    # attrdict creates a copy of dct.
    ad.a = 2
    assert ad.a == 2 and dct["a"] == 1
    ad.c.d = 4
    assert ad.c.d == 4 and dct["c"]["d"] == 3
