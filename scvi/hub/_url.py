import re

import requests


def validate_url(url: str, raise_error: bool = False) -> bool:
    """Validates a URL.

    Source: https://stackoverflow.com/questions/7160737/how-to-validate-a-url-in-python-malformed-or-not
    """
    regex = re.compile(
        r"^(?:http|ftp)s?://"  # http:// or https://
        r"(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|"  # domain...
        r"localhost|"  # localhost...
        r"\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})"  # ...or ip
        r"(?::\d+)?"  # optional port
        r"(?:/?|[/?]\S+)$",
        re.IGNORECASE,
    )
    if re.match(regex, url) is None:
        if raise_error:
            raise ValueError(f"Invalid URL format: {url}")
        return False

    try:
        response = requests.get(url)
        valid = response.status_code == 200
    except requests.ConnectionError:
        valid = False

    if not valid and raise_error:
        raise ValueError(f"Invalid URL: {url}")

    return valid
