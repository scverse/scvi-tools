import logging
import os
import urllib

import numpy as np

from scvi.utils import track

logger = logging.getLogger(__name__)


def _download(url: str | None, save_path: str, filename: str):
    """Writes data from url to file."""
    download_link = os.path.join(save_path, filename)
    if os.path.exists(download_link):
        logger.info(f"File {download_link} already downloaded")
        return
    elif url is None:
        logger.info(f"No backup URL provided for missing file {download_link}")
        return
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
