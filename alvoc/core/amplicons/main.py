from pathlib import Path
import pandas as pd
from alvoc.core.amplicons.visualize import plot_depths, plot_depths_gc
from alvoc.core.utils import create_dir, precompute
import pysam


def calculate_amplicon_metrics(
    virus: str,
    samples: Path,
    inserts: Path,
    outdir: Path,
    max_depth: int = 50000,
) -> pd.DataFrame:
    """
    Calculates the amplicon coverage and GC content for samples listed in a BAM file or a samplesheet.

    Args:
        virus: Name of the virus.
        samples: Path to a BAM file or CSV file listing sample BAM paths and labels.
        inserts: Path to the inserts file (BED format).
        outdir: Path to the output directory.
        max_depth: Maximum depth to be considered for calculations.

    Returns:
        A Pandas DataFrame containing amplicon identifiers, read depths, and GC content for each sample.
    """
    # Create output directory
    out = create_dir(outdir=outdir)

    # Run precompute to get the sequence
    _, seq = precompute(virus, out)

    # Parse BED file to extract amplicon information
    with open(inserts, "r") as bed:
        inserts_list = [
            line.strip().split("\t")[:4]
            for line in bed.readlines()
            if len(line.strip().split("\t")) >= 4
        ]

    all_results = []

    if samples.suffix == ".bam":
        # Process single BAM file
        sample_name = samples.stem
        coverage = find_depths_in_bam(samples, inserts_list, max_depth=max_depth)
        gc_content = calculate_gc_depth(seq, inserts_list)

        # Combine results for this sample
        all_results.extend(
            [
                {
                    "sample": sample_name,
                    "amplicon_id": amplicon,
                    "depth": coverage[amplicon],
                    "gc_content": gc_content[amplicon],
                }
                for amplicon in coverage
            ]
        )

    elif samples.suffix == ".csv":
        # Process multiple samples from a CSV file
        sample_df = pd.read_csv(samples)
        for _, row in sample_df.iterrows():
            bam_path = Path(row["bam"])
            sample_name = row["sample"]
            if bam_path.suffix == ".bam":
                coverage = find_depths_in_bam(
                    bam_path, inserts_list, max_depth=max_depth
                )
                gc_content = calculate_gc_depth(seq, inserts_list)

                # Combine results for this sample
                all_results.extend(
                    [
                        {
                            "sample": sample_name,
                            "amplicon_id": amplicon,
                            "depth": coverage[amplicon],
                            "gc_content": gc_content[amplicon],
                        }
                        for amplicon in coverage
                    ]
                )
    results_df = pd.DataFrame(all_results)

    # Save results to CSV
    if not results_df.empty:
        results_df.to_csv(out / "amplicon_metrics.csv", index=False)
        plot_depths(results_df, inserts_list, out)
        plot_depths_gc(results_df, out)

    return results_df


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


def calculate_gc_depth(sequence: str, inserts: list[list]) -> dict:
    """
    Calculates the GC content for each amplicon region.

    Args:
        sequence: The full genome sequence as a string.
        inserts: List of Lists containing information about the regions (amplicons) of interest.

    Returns:
        A dictionary mapping amplicon identifiers to their corresponding GC content.
    """
    gc_content = {}
    for insert in inserts:
        start, end, amplicon_id = int(insert[1]), int(insert[2]), insert[3]
        region_seq = sequence[start:end]
        gc_count = region_seq.count("G") + region_seq.count("C")
        gc_content[amplicon_id] = (
            gc_count / len(region_seq) if len(region_seq) > 0 else 0
        )

    return gc_content
