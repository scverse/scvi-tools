import logging
import os
from typing import Optional

import scanpy as sc

logger = logging.getLogger(__name__)


def _download(url: Optional[str], save_path: str, filename: str):
    """Writes data from url to file."""
    path_filename = os.path.join(save_path, filename)
    if os.path.exists(path_filename):
        logger.info(f"File {path_filename} already downloaded")
        return
    elif url is None:
        logger.info(f"No backup URL provided for missing file {path_filename}")
        return
    sc.readwrite._check_datafile_present_and_download(path_filename, url)
