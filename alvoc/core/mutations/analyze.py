from pathlib import Path
from alvoc.core.mutations.utils import (
    mut_idx,
    write_csv,
    snv_name,
    mut_in_col,
    parse_mutation,
    print_mut_results,
)
from alvoc.core.mutations.visualize import plot_mutations
import pysam


def find_mutants(
    file_path: str,
    mutations_path: str,
    min_depth: int, 
    mut_lins: dict,
    genes: dict,
    seq: str,
    outdir: Path,
):
    """Find mutations in sequencing data, either from BAM files or a sample list. Uses a dictionary of mutation lineages provided as a parameter.

    Args:
        file_path (str): Path to the file containing sample information or BAM file.
        mutations_path (str): Path to the file containing mutations or mutation identifier.
        min_depth (int): Minimum depth for mutation analysis.
        save_img (bool): Whether to save a plot image.
        csv (bool): Whether to generate a CSV file.
        mut_lins (dict): Dictionary containing mutation lineages and their occurrences.

    Returns:
        None: The function directly modifies files and outputs results.
    """
    sample_results, sample_names = [], [] 
     # Function to adapt mut_idx for sorting
    def mut_idx_adapter(mut):
        return mut_idx(mut, genes, seq)

    # Determine if mutations_path is a known lineage or a file with mutations
    if mutations_path in mut_lins:
        print("Searcing for {} mutations".format(mutations_path))
        mutations = [
            mut
            for mut in mut_lins
            if mut_lins[mut][mutations_path] > 0 and mut_idx(mut, genes, seq) != -1
        ]
        mutations.sort(key=mut_idx_adapter)
    else:
        with open(mutations_path, "r") as file:
            mutations = [mut.strip() for mut in file.read().split("\n") if mut.strip()]

    # Handle BAM files or sample lists
    if file_path.endswith(".bam"):
        sample_results.append(find_mutants_in_bam(file_path, mutations, genes, seq))
        sample_names.append("")
    else:
        with open(file_path, "r") as file:
            samples = [line.split("\t") for line in file.read().split("\n") if line]
        for sample in samples:
            if sample[0].endswith(".bam"):
                sample_results.append(find_mutants_in_bam(sample[0], mutations, genes, seq))
                sample_names.append(sample[1])
                print_mut_results(sample_results[-1], min_depth)

    mutants_name = mutations_path.rsplit(".", 1)[0]
    write_csv(sample_results, sample_names, min_depth, mutants_name, outdir)
    plot_mutations(sample_results, sample_names, min_depth, mutants_name, outdir)


def find_mutants_in_bam(bam_path, mutations, genes, seq):
    """Identify and quantify mutations from a BAM file.

    Args:
        bam_path (str): Path to the BAM file.
        mutations (list): A list of mutations to look for in the BAM file.

    Returns:
        dict: A dictionary where keys are mutations and values are the maximal frequency
              of each mutation and its count of occurrences and non-occurrences.
    """
    mut_results = {}

    with pysam.Samfile(bam_path, "rb") as samfile:
        parsed_muts = {mut: parse_mutation(mut, genes, seq) for mut in mutations}
        mut_results = {
            mut: {snv_name(m): [0, 0] for m in parsed_muts[mut]} for mut in parsed_muts
        }

        # Iterate over each pileup column in the BAM file
        for pileupcolumn in samfile.pileup(stepper="nofilter"):
            pos = pileupcolumn.reference_start + 1
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


def evaluate_mutation_frequencies(mut_results):
    """Evaluate the frequency of each mutation in the results.

    Args:
        mut_results (dict): A dictionary containing counts of mutations.

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
