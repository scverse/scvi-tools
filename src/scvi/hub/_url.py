import re

import requests

from scvi.utils import dependencies


def validate_url(url: str, error_format: bool = False, error_response: bool = False) -> bool:
    """Validates a URL.

    Source: https://stackoverflow.com/questions/7160737/how-to-validate-a-url-in-python-malformed-or-not
    """
    regex = re.compile(
        r"^(?:http|ftp)s?://"  # http:// or https://
        r"(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|"
        r"localhost|"  # localhost...
        r"\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})"  # ...or ip
        r"(?::\d+)?"  # optional port
        r"(?:/?|[/?]\S+)$",
        re.IGNORECASE,
    )
    if re.match(regex, url) is None:
        if error_format:
            raise ValueError(f"Invalid URL format: {url}")
        return False

    try:
        response = requests.get(url)
        valid = response.status_code == 200
    except requests.ConnectionError:
        valid = False

    if not valid and error_response:
        raise ValueError(f"Invalid URL: {url}")

    return valid


@dependencies("selenium")
def validate_colab_notebook(colab_url: str) -> bool:
    from selenium import webdriver
    from selenium.webdriver.chrome.options import Options

    timeout = 15
    options = Options()
    options.add_argument("--headless")
    options.add_argument("--disable-gpu")
    options.add_argument("--no-sandbox")

    driver = webdriver.Chrome(options=options)

    driver.get(colab_url)

    # Wait for content to load
    driver.implicitly_wait(timeout)

    page_source = driver.page_source

    if "Notebook not found" in page_source:
        print(f"❌ Notebook not found: {colab_url}")
        driver.quit()
        return False
    elif "scvi-tools" in page_source:
        print(f"✅ Valid notebook: {colab_url}")
        driver.quit()
        return True
    else:
        print(f"⚠️ Unknown state: {colab_url}")
        driver.quit()
        return False
