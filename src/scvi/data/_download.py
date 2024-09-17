import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)


def _download_scanpy_style(url: str, path: Path):
    try:
        import ipywidgets  # noqa: F401
        from tqdm.auto import tqdm
    except ImportError:
        from tqdm import tqdm

    from urllib.error import URLError
    from urllib.request import Request, urlopen

    blocksize = 1024 * 8
    blocknum = 0

    try:
        req = Request(url, headers={"User-agent": "scanpy-user"})

        try:
            open_url = urlopen(req)
        except URLError:
            logger.warning(
                "Failed to open the url with default certificates, trying with certifi."
            )

            from ssl import create_default_context

            from certifi import where

            open_url = urlopen(req, context=create_default_context(cafile=where()))

        with open_url as resp:
            total = resp.info().get("content-length", None)
            with (
                tqdm(
                    unit="B",
                    unit_scale=True,
                    miniters=1,
                    unit_divisor=1024,
                    total=total if total is None else int(total),
                ) as t,
                path.open("wb") as f,
            ):
                block = resp.read(blocksize)
                while block:
                    f.write(block)
                    blocknum += 1
                    t.update(len(block))
                    block = resp.read(blocksize)

    except (KeyboardInterrupt, Exception):
        # Make sure file doesn’t exist half-downloaded
        if path.is_file():
            path.unlink()
        raise


def _check_datafile_present_and_download(path, backup_url=None):
    """Check whether the file is present, otherwise download."""
    path = Path(path)
    if path.is_file():
        return True
    if backup_url is None:
        return False
    logger.info(
        f"try downloading from url\n{backup_url}\n"
        "... this may take a while but only happens once"
    )
    if not path.parent.is_dir():
        logger.info(f"creating directory {path.parent}/ for saving data")
        path.parent.mkdir(parents=True)

    _download_scanpy_style(backup_url, path)
    return True


def _download(url: str | None, save_path: str, filename: str):
    """Writes data from url to file."""
    path_filename = os.path.join(save_path, filename)
    if os.path.exists(path_filename):
        logger.info(f"File {path_filename} already downloaded")
        return
    elif url is None:
        logger.info(f"No backup URL provided for missing file {path_filename}")
        return
    _check_datafile_present_and_download(path_filename, url)
