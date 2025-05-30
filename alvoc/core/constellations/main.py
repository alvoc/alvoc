from alvoc.core.constellations.core import (
    DataSource,
    calculate_profiles,
    create_constellations,
)
from pathlib import Path


def make_constellations(
    source: DataSource, source_path: str, outdir: Path, threshold: float
):
    raw = source.fetch(source_path)
    profiles = calculate_profiles(source, raw, threshold)
    create_constellations(profiles, outdir)
    print("Done.")
