import logging
import os
import time
import urllib

import numpy as np

from scvi.utils import track

logger = logging.getLogger(__name__)

# transient failures worth a fresh attempt: dropped/reset connections partway through a
# download, as opposed to e.g. a 404 which will never succeed on retry
_TRANSIENT_DOWNLOAD_ERRORS = (
    ConnectionError,
    TimeoutError,
    urllib.error.URLError,
)
_MAX_DOWNLOAD_ATTEMPTS = 3
_RETRY_BACKOFF_SECONDS = 5


def _download(url: str | None, save_path: str, filename: str):
    """Writes data from url to file."""
    download_link = os.path.join(save_path, filename)
    if os.path.exists(download_link):
        logger.info(f"File {download_link} already downloaded")
        return
    elif url is None:
        logger.info(f"No backup URL provided for missing file {download_link}")
        return

    def read_iter(file, block_size=1000):
        """Iterates through file.

        Given a file 'file', returns an iterator that returns bytes of
        size 'blocksize' from the file, using read().
        """
        while True:
            block = file.read(block_size)
            if not block:
                break
            yield block

    # Create the path to save the data
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    block_size = 1000

    for attempt in range(1, _MAX_DOWNLOAD_ATTEMPTS + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Magic Browser"})
            r = urllib.request.urlopen(req)
            if (r.getheader("Content-Length") is None) and (
                r.getheader("Content-Type") != "text/tab-separated-values"
            ):
                raise FileNotFoundError(
                    f"Found file with no content at {url}. "
                    "This is possibly a directory rather than a file path."
                )

            logger.info(f"Downloading file at {download_link}")
            if r.getheader("Content-Length") is not None:
                filesize = int(r.getheader("Content-Length"))
                filesize = np.rint(filesize / block_size)
                with open(download_link, "wb") as f:
                    iterator = read_iter(r, block_size=block_size)
                    for data in track(
                        iterator, style="tqdm", total=filesize, description="Downloading..."
                    ):
                        f.write(data)
            else:
                urllib.request.urlretrieve(url, download_link)
                print(f"File downloaded successfully and saved as {download_link}")
            return
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raise FileNotFoundError(f"Could not find file at {url}") from exc
            if attempt == _MAX_DOWNLOAD_ATTEMPTS or exc.code < 500:
                raise
            _wait_before_retry(url, attempt, exc)
        except _TRANSIENT_DOWNLOAD_ERRORS as exc:
            if os.path.exists(download_link):
                os.remove(download_link)
            if attempt == _MAX_DOWNLOAD_ATTEMPTS:
                raise
            _wait_before_retry(url, attempt, exc)


def _wait_before_retry(url: str, attempt: int, exc: Exception) -> None:
    logger.warning(
        f"Download of {url} was interrupted (attempt {attempt}/"
        f"{_MAX_DOWNLOAD_ATTEMPTS}): {exc!r}. Retrying in {_RETRY_BACKOFF_SECONDS}s."
    )
    time.sleep(_RETRY_BACKOFF_SECONDS)


def _pooch_retrieve_with_retries(**retrieve_kwargs):
    """Calls ``pooch.retrieve`` with retries for transient network failures.

    ``retry_if_failed`` is only accepted by ``pooch.create``/``Pooch.fetch``, not by the
    one-off ``pooch.retrieve`` helper used here: passing it to ``pooch.retrieve`` (or to
    ``HTTPDownloader``) is silently forwarded as far as ``requests.get``, which raises
    ``TypeError`` for the unknown keyword. So retries are implemented here instead, by
    re-issuing the whole ``pooch.retrieve`` call, which works across pooch versions.
    """
    import pooch
    import requests

    for attempt in range(1, _MAX_DOWNLOAD_ATTEMPTS + 1):
        try:
            return pooch.retrieve(**retrieve_kwargs)
        except (requests.exceptions.RequestException, *_TRANSIENT_DOWNLOAD_ERRORS) as exc:
            if attempt == _MAX_DOWNLOAD_ATTEMPTS:
                raise
            _wait_before_retry(retrieve_kwargs.get("url"), attempt, exc)
