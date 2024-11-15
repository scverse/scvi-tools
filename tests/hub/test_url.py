from scvi.hub._url import validate_url


def test_validate_url():
    valid = "https://scvi-tools.org/"
    assert validate_url(valid)

    invalid = "scvi-tools.org"
    assert not validate_url(invalid)
