from pathlib import Path
from typing import Callable
import pysam
from alvoc.core.amplicons.visualize import plot_depths, plot_depths_gc
import pandas as pd


def amplicon_coverage(file_path: Path, inserts: Path, sequence: str, outdir: Path):
    """
    Determines and plots amplicon coverage for samples listed in a file or a single BAM file.

    Args:
        file_path: Path to the file listing samples or a single BAM file.
                        If a file listing samples, it should contain one path to a BAM file per line.
        inserts: Path to the CSV file detailing the regions (amplicons) to evaluate.
                        The CSV file should have columns: 'chrom', 'chromStart', 'chromEnd'.
    """
    process_samples(file_path, inserts, plot_depths, sequence, outdir)


def gc_depth(file_path: Path, inserts: Path, sequence: str, outdir: Path):
    """
    Determines and plots the GC depth correlation for samples listed in a file or a single BAM file.

    Args:
        file_path: Path to the file listing samples or a single BAM file.
                        If a file listing samples, it should contain one path to a BAM file per line.
        inserts: Path to the CSV file detailing the regions (amplicons) to evaluate.
                        The CSV file should have columns: 'chrom', 'chromStart', 'chromEnd'.
    """
    process_samples(file_path, inserts, plot_depths_gc, sequence, outdir)


def process_samples(
    file_path: Path, inserts: Path, plot_function: Callable, sequence: str, outdir: Path
):
    """
    Processes a file containing sample paths and names, extracts depths from BAM files, and visualizes them.

    Args:
        file_path: Path to the file listing samples or a single BAM file.
                        If a file listing samples, it should contain one path to a BAM file per line.
        inserts: Path to the CSV file detailing the regions (amplicons) to evaluate.
                        The CSV file should have columns: 'chrom', 'chromStart', 'chromEnd'.

    Examples:
        Example of a text file listing sample BAM files:
        ```
        /path/to/sample1.bam Sample 1
        /path/to/sample2.bam Sample 2
        ```

        Example of a CSV file detailing regions:
        ```
        chromosome,start,end
        chr1,10000,10500
        chr2,20000,20500
        ```
    """
    bed_df = pd.read_csv(inserts)

    gc_contents = []
    for _, row in bed_df.iterrows():
        section = sequence[int(row["chromStart"]) : int(row["chromEnd"])]
        gc_content = calculate_gc_content(section)
        gc_contents.append(gc_content)

    bed_df["gcContent"] = gc_contents

    final_inserts = bed_df.values.tolist()
    sample_results, sample_names = [], []
    if file_path.suffix == ".bam":
        sample_results.append(find_depths_in_bam(file_path, final_inserts))
        sample_names.append("")
    else:
        with file_path.open("r") as f:
            samples = [
                line.strip().split("\t") for line in f.readlines() if line.strip()
            ]
        for bam, sample in samples:
            bam_path = Path(bam)
            if bam_path.suffix == ".bam":
                sample_results.append(find_depths_in_bam(bam_path, final_inserts))
                sample_names.append(sample)
    plot_function(sample_results, sample_names, inserts, outdir)


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate the GC content of a given sequence.

    Args:
        sequence: DNA sequence.

    Returns:
        GC content as a fraction.
    """
    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence) if len(sequence) > 0 else 0


def find_depths_in_bam(
    bam_path: Path, inserts: list[list], max_depth: int = 50000
) -> dict:
    """
    Reads a BAM file and computes the depth of reads at positions defined by `inserts`.

    Args:
    bam_path: Path to the BAM file.
    inserts: List of Lists containing information about the regions (amplicons) of interest.
    max_depth: Maximum depth to be considered to prevent memory overflow.

    Returns:
    A dictionary mapping amplicon identifiers to their corresponding read depth.
    """
    amplified = {}
    with pysam.AlignmentFile(bam_path.as_posix(), "rb") as samfile:
        amp_mids = {int((int(i[1]) + int(i[2])) / 2): i[3] for i in inserts}
        amplified = {i[3]: 0 for i in inserts}
        for pileupcolumn in samfile.pileup(max_depth=max_depth):
            pos = pileupcolumn.reference_pos
            if pos in amp_mids:
                depth = pileupcolumn.get_num_aligned()
                amplified[amp_mids[pos]] = depth

    return amplified
