from pathlib import Path
from alvoc.core.utils import logging

logger = logging.get_logger()


def create_dir(outdir: Path):
    outdir_path = Path(outdir)
    if not outdir_path.is_dir():
        logger.info(f"Creating directory at {outdir_path}")
        outdir_path.mkdir(parents=True, exist_ok=True)
    return outdir_path
