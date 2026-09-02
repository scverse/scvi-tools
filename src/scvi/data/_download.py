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
        req = urllib.request.Request(url, headers={"User-Agent": "Magic Browser"})
        try:
            r = urllib.request.urlopen(req)
            if (r.getheader("Content-Length") is None) and (
                r.getheader("Content-Type") != "text/tab-separated-values"
            ):
                raise FileNotFoundError(
                    f"Found file with no content at {url}. "
                    "This is possibly a directory rather than a file path."
                )
        except urllib.error.HTTPError as exc:
            if exc.code == "404":
                raise FileNotFoundError(f"Could not find file at {url}") from exc
            raise exc

        logger.info(f"Downloading file at {download_link}")
        try:
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
        except _TRANSIENT_DOWNLOAD_ERRORS:
            if os.path.exists(download_link):
                os.remove(download_link)
            if attempt == _MAX_DOWNLOAD_ATTEMPTS:
                raise
            logger.warning(
                f"Download of {url} was interrupted (attempt {attempt}/"
                f"{_MAX_DOWNLOAD_ATTEMPTS}), retrying in {_RETRY_BACKOFF_SECONDS}s."
            )
            time.sleep(_RETRY_BACKOFF_SECONDS)
