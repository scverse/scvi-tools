import logging
import os
import urllib

import numpy as np

from scvi._utils import track

logger = logging.getLogger(__name__)


def _download(url: str, save_path: str, filename: str):
    """Writes data from url to file."""
    if os.path.exists(os.path.join(save_path, filename)):
        logger.info("File %s already downloaded" % (os.path.join(save_path, filename)))
        return
    req = urllib.request.Request(url, headers={"User-Agent": "Magic Browser"})
    r = urllib.request.urlopen(req)
    logger.info("Downloading file at %s" % os.path.join(save_path, filename))

    def read_iter(file, block_size=1000):
        """
        Iterates through file.

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

    filesize = int(r.getheader("Content-Length"))
    filesize = np.rint(filesize / block_size)
    with open(os.path.join(save_path, filename), "wb") as f:
        iterator = read_iter(r, block_size=block_size)
        for data in track(
            iterator, style="tqdm", total=filesize, description="Downloading..."
        ):
            f.write(data)
