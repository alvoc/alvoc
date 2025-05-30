from abc import ABC, abstractmethod
from pathlib import Path
import json
from collections import defaultdict
from typing import Any, Iterator, List, Tuple


class DataSource(ABC):
    """
    Abstract base class for phylogenetic data sources.
    Each source must implement `fetch` to load raw data and
    `records` to yield (clade_name, mutations_list) tuples.
    """

    @abstractmethod
    def fetch(self, source: str) -> object:
        """
        Download or load the raw data from a URL or filepath.
        Returns a source-specific raw data object.
        """
        pass

    @abstractmethod
    def records(self, raw_data: Any) -> Iterator[Tuple[str, List[str]]]:
        """
        Iterate over all entries, yielding:
          clade_name (str), mutations (list of str)
        """
        ...


# Generic pipeline

def calculate_profiles(source: DataSource, raw_data: object, threshold: float):
    """
    Generic function to compute defining mutations per clade from any DataSource.
    Returns a dict: clade_name → set of mutations.
    """
    # count sequences/nodes and mutation occurrences
    clade_counts = defaultdict(int)
    mut_counts = defaultdict(lambda: defaultdict(int))

    for clade, muts in source.records(raw_data):
        clade_counts[clade] += 1
        for m in muts:
            mut_counts[clade][m] += 1

    # filter by threshold
    profiles = {}
    for clade, counts in mut_counts.items():
        total = clade_counts[clade]
        kept = {m for m, c in counts.items() if (c / total) >= threshold}
        if kept:
            profiles[clade] = kept

    return profiles


def create_constellations(profiles: dict, output_dir: Path):
    """Write out 'constellations.json' and 'lineages.txt'."""
    const = {}
    for clade, muts in profiles.items():
        const[clade] = {
            "lineage": clade,
            "label": f"{clade}-like",
            "description": f"{clade} lineage defining mutations",
            "sources": [],
            "tags": [clade],
            "sites": list(muts),
            "note": "Unique mutations for sublineage",
        }
    # JSON
    out_json = output_dir / "constellations.json"
    with open(out_json, "w") as f:
        json.dump(const, f, indent=4)
    # TXT
    out_txt = output_dir / "lineages.txt"
    with open(out_txt, "w") as f:
        for clade in sorted(const):
            f.write(clade + "\n")
    return const


def make_constellations(
    source: DataSource, source_path: str, outdir: Path, threshold: float
):
    raw = source.fetch(source_path)
    profiles = calculate_profiles(source, raw, threshold)
    create_constellations(profiles, outdir)
    print("Done.")
