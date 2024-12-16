from pathlib import Path

import pandas as pd
import pysam

from alvoc.core.variants.mutations.helpers import mut_in_col
from alvoc.core.variants.mutations.visualize import plot_mutations
from alvoc.core.variants.prepare import parse_lineages
from alvoc.core.utils.parse import mut_idx, parse_mutation, snv_name
from alvoc.core.utils import create_dir, logging, precompute

logger = logging.get_logger()


def find_mutants(
    virus: str,
    samples: Path,
    constellations: Path,
    mutations_path: Path,
    outdir: Path,
    min_depth: int = 10,
):
    """Find mutations in sequencing data, either from BAM files or a sample list. Uses a dictionary of mutation lineages provided as a parameter.

    Args:
        virus: Taxonomic ID of the virus, or path to the GenBank file
        samples: Path to a BAM file or CSV file listing samples.
        constellations: Path to a JSON file containing mutation lineage constellations.
        mutations_path: Path to the file containing mutations or mutation identifier.
        min_depth: Minimum depth for mutation analysis.
        mut_lins: Dictionary containing mutation lineages and their occurrences.
        outdir : Output directory for results and intermediate data. Defaults to the current directory.

    Returns:
        None: The function directly modifies files and outputs results.
    """
    # Create or find directory
    out = create_dir(outdir)

    # Extract the genome sequence and gene coordinates for the target virus
    genes, seq = precompute(virus, out)

    # Convert lineage data to mutation-centric format
    mut_lins = parse_lineages(constellations)

    results_df = pd.DataFrame()

    # Function to adapt mut_idx for sorting
    def mut_idx_adapter(mut):
        return mut_idx(mut, genes, seq)

    # Determine if mutations_path is a known lineage or a file with mutations
    if mutations_path.name in mut_lins:
        print(f"Searching for {mutations_path.name} mutations")
        mutations = [
            mut
            for mut in mut_lins
            if mut_lins[mut][mutations_path.name] > 0 and mut_idx(mut, genes, seq) != -1
        ]
        mutations.sort(key=mut_idx_adapter)
    else:
        with mutations_path.open("r") as file:
            mutations = [mut.strip() for mut in file.read().split("\n") if mut.strip()]

    if samples.suffix == ".bam":
        results = find_mutants_in_bam(
            bam_path=samples,
            mutations=mutations,
            genes=genes,
            seq=seq,
        )
        mr_df = pd.DataFrame.from_dict(
            results,
            orient="index",
            columns=["mutation_count", "non_mutation_count"],
        ).reset_index(names=["mutants"])
        mr_df.insert(1, "sample", samples.stem)
        results_df = pd.concat([results_df, mr_df], ignore_index=True)
    else:
        sample_df = pd.read_csv(samples)
        for _, row in sample_df.iterrows():
            bam_path = Path(row["bam"])
            sample_label = row["sample"]
            if bam_path.suffix == ".bam":
                results = find_mutants_in_bam(
                    bam_path=bam_path,
                    mutations=mutations,
                    genes=genes,
                    seq=seq,
                )
                mr_df = pd.DataFrame.from_dict(
                    results,
                    orient="index",
                    columns=["mutation_count", "non_mutation_count"],
                ).reset_index(names=["mutants"])
                mr_df.insert(1, "sample", sample_label)
                results_df = pd.concat([results_df, mr_df], ignore_index=True)

    if not results_df.empty:
        results_df.to_csv(out / "mutations_melted.csv", index=False)
        plot_mutations(results_df, min_depth, out)


def find_mutants_in_bam(bam_path: Path, mutations, genes, seq):
    """Identify and quantify mutations from a BAM file.

        Args:
            bam_path (Path): Path to the BAM file.
            mutations (list): A list of mutations to look for in the BAM file.

        Returns:
            dict: A dictionary where keys are mutation names, and values are lists containing:
                - The count of occurrences (`mutation_count`) of the single nucleotide variation (SNV) with the highest observed frequency.
                - The count of non-occurrences (`non_mutation_count`) corresponding to that SNV.
    .
    """
    mut_results = {}

    with pysam.Samfile(bam_path.as_posix(), "rb") as samfile:
        parsed_muts = {mut: parse_mutation(mut, genes, seq) for mut in mutations}
        mut_results = {
            mut: {snv_name(m): [0, 0] for m in parsed_muts[mut]} for mut in parsed_muts
        }

        # Iterate over each pileup column in the BAM file
        for pileupcolumn in samfile.pileup(stepper="nofilter"):
            pos = pileupcolumn.reference_pos + 1
            update_mutation_results(pileupcolumn, parsed_muts, mut_results, pos)

    output = evaluate_mutation_frequencies(mut_results)
    return output


def update_mutation_results(pileupcolumn, parsed_muts, mut_results, pos):
    """Update mutation results based on pileup column data.

    Args:
        pileupcolumn (PileupColumn): Pileup column object from a BAM file.
        parsed_muts (dict): A dictionary containing parsed mutations.
        mut_results (dict): A dictionary to store results of mutation counts.
        pos (int): Current position in the BAM file being examined.

    Returns:
        None: Modifies `mut_results` in place.
    """
    for mut, snvs in parsed_muts.items():
        for snv in snvs:
            if pos == snv[1]:  # Check if position matches the mutation position
                muts, not_muts = mut_in_col(pileupcolumn, snv[2])
                mut_results[mut][snv_name(snv)] = [muts, not_muts]


def evaluate_mutation_frequencies(mut_results: dict):
    """Evaluate the frequency of each mutation in the results.

    Args:
        mut_results: A dictionary containing counts of mutations.

    Returns:
        dict: A dictionary with each mutation and its highest observed frequency.
    """
    for mut, results in mut_results.items():
        max_freq = -1
        max_muts = [0, 0]
        for result in results.values():
            muts, not_muts = result
            total = muts + not_muts
            freq = muts / total if total > 0 else 0
            if freq > max_freq:
                max_freq = freq
                max_muts = [muts, not_muts]
        mut_results[mut] = max_muts
    return mut_results
