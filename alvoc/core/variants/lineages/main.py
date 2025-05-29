from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd


from alvoc.core.utils import create_dir, logging, precompute
from alvoc.core.variants.lineages.regression import do_regression, do_regression_linear
from alvoc.core.variants.mutations import find_mutants_in_bam
from alvoc.core.variants.prepare import parse_lineages
from alvoc.core.variants.lineages.visualize import (
    plot_lineages,
)
from alvoc.core.utils.parse import parse_mutations

logger = logging.get_logger()


def find_lineages(
    virus: str,
    samples: Path,
    constellations: Path,
    outdir: Path,
    white_list: Path | None = None,
    black_list: Path | None = None,
    min_depth: int = 40,
    unique: bool = False,
    l2: bool = False,
):
    """
    Processes either a single BAM file or a samplesheet to find and analyze lineages based on the given parameters.

    Args:
        virus: Taxonomic ID of the virus, or path to the GenBank file
        samples: Path to a BAM file or CSV file listing samples.
        constellations: Path to a JSON file containing mutation lineage constellations.
        outdir: Output directory for results and intermediate data. Defaults to the current directory.
        white_list: Path to a TXT file containing lineages to include.
        black_list: Path to a TXT file containing lineages to exclude.
        min_depth: Minimum depth for a mutation to be considered.
        unique: Whether to consider unique mutations only.
        l2: Whether to use a secondary method for regression analysis.

    Returns:
        None: Outputs csv + plots.
    """
    # Create or find directory
    out = create_dir(outdir)

    # Extract the genome sequence and gene coordinates for the target virus
    genes, seq = precompute(virus, out)

    # Convert lineage data to mutation-centric format
    mut_lins = parse_lineages(constellations)

    results_df = pd.DataFrame()

    if samples.suffix == ".bam":
        # Process single BAM file
        results_df = process_sample(
            sample=samples,
            sample_name=samples.stem,
            results_df=results_df,
            mut_lins=mut_lins,
            genes=genes,
            seq=seq,
            white_list=white_list,
            black_list=black_list,
            min_depth=min_depth,
            unique=unique,
            l2=l2,
        )
    else:
        # Process multiple samples from a CSV file
        sample_df = pd.read_csv(samples)
        for _, row in sample_df.iterrows():
            bam_path = Path(row["bam"])
            sample_label = row["sample"]
            if bam_path.suffix == ".bam":
                results_df = process_sample(
                    sample=bam_path,
                    sample_name=sample_label,
                    results_df=results_df,
                    mut_lins=mut_lins,
                    genes=genes,
                    seq=seq,
                    white_list=white_list,
                    black_list=black_list,
                    min_depth=min_depth,
                    unique=unique,
                    l2=l2,
                )

    if not results_df.empty:
        results_df.to_csv(out / "lineage_abundance_melted.csv", index=False)
        pivoted_df = results_df.pivot(
            index="lineage", columns="sample", values="abundance"
        )
        pivoted_df.to_csv(out / "lineage_abundance_pivoted.csv")
        plot_lineages(results_df, out, bool(white_list))


def process_sample(
    sample: Path,
    sample_name: str,
    results_df: pd.DataFrame,
    mut_lins: dict,
    genes: dict,
    seq: str,
    white_list: Path | None = None,
    black_list: Path | None = None,
    min_depth: int = 40,
    unique: bool = False,
    l2: bool = False,
) -> pd.DataFrame:
    """Quantify lineages for a single BAM file."""
    # Load white list lineages if provided
    included_lineages = []
    if white_list:
        with open(white_list, "r") as f:
            included_lineages = f.read().splitlines()

    # Load black listed lineages if provided
    excluded_lineages = []
    if black_list:
        with open(black_list, "r") as f:
            excluded_lineages = f.read().splitlines()

    result = quantify_lineages(
        sample=sample,
        mut_lins=mut_lins,
        genes=genes,
        seq=seq,
        white_list=included_lineages,
        black_list=excluded_lineages,
        min_depth=min_depth,
        unique=unique,
        l2=l2,
    )
    if result is None:
        logger.error(
            f"No coverage or analysis couldn't be performed for {sample_name}."
        )
        return results_df

    sr, _, _, _ = result
    sr_df = pd.DataFrame(list(sr.items()), columns=["lineage", "abundance"])
    sr_df.insert(0, "sample", sample_name)
    return pd.concat([results_df, sr_df], ignore_index=True)


def quantify_lineages(
    sample: Path,
    mut_lins: dict,
    genes: dict,
    seq: str,
    black_list: list[str] = [],
    white_list: list[str] = [],
    min_depth: int = 40,
    unique: bool = False,
    l2: bool = False,
) -> Union[None, tuple[dict[str, float], np.ndarray, np.ndarray, list[str]]]:
    """
    Identify and estimate abundance of lineages in a BAM file based on predefined mutations.

    Args:
        samples: Path to a BAM file or CSV file listing samples.
        mut_lins: Dictionary containing mutation lineages and their occurrences.
        genes: Dictionary mapping gene names to their start and end positions in the sequence.
        seq: Complete nucleotide sequence as a string.
        outdir: Output directory for results and intermediate data. Defaults to the current directory.
        white_list: Path to a TXT file containing lineages to include.
        black_list: Path to a TXT file containing lineages to exclude.
        min_depth: Minimum depth for a mutation to be considered.
        unique: Whether to consider unique mutations only.
        l2: Whether to use a secondary method for regression analysis.

    Returns:
        Sample results, optionally with additional data structures.
    """
    # Make all mutation and lineage names uppercase for case-insensitive matching
    mut_lins = {m.upper(): {lin.upper(): v for lin, v in lins.items()} for m, lins in mut_lins.items()}
    black_list = [b.upper() for b in black_list]
    white_list = [w.upper() for w in white_list]

    logger.info("Identifying target lineages")
    aa_mutations = [m for m in mut_lins.keys() if m not in black_list]
    if unique:
        aa_mutations = [
            mut
            for mut in aa_mutations
            if sum(mut_lins[mut][lin] for lin in white_list) == 1
        ]

    logger.info("Converting to nucleotides")
    mutations = parse_mutations(aa_mutations, genes, seq)

    logger.info("Finding mutants")
    mut_results = find_mutants_in_bam(
        bam_path=sample, mutations=mutations, genes=genes, seq=seq
    )

    logger.info("Filtering out mutations below min_depth")
    covered_muts = [m for m in mutations if sum(mut_results[m]) >= min_depth]
    if not covered_muts:
        logger.info("No coverage")
        return None

    covered_lineages = {
        lin
        for m in covered_muts
        for lin in mut_lins[m]
        if mut_results[m][0] > 0 and mut_lins[m][lin] > 0.5
    }

    if white_list:
        covered_lineages = [lin for lin in white_list if lin in covered_lineages]
    else:
        covered_lineages = list(covered_lineages)

    Y = np.array([mut_results[m][0] / sum(mut_results[m]) for m in covered_muts])
    lmps = [[round(mut_lins[m][lin]) for m in covered_muts] for lin in covered_lineages]

    # Merge indistinguishable lineages to accurately reflect diversity
    merged_lmps = []
    merged_lins = []
    for i in range(len(covered_lineages)):
        lmp = lmps[i]
        lin = covered_lineages[i]
        if lmp in merged_lmps:
            lmp_idx = merged_lmps.index(lmp)
            merged_lins[lmp_idx] += " or " + lin
        else:
            merged_lmps.append(lmp)
            merged_lins.append(lin)

    if l2:
        X, reg = do_regression(lmps, Y)
    else:
        X, reg, mut_diffs = do_regression_linear(lmps, Y, covered_muts)
    
    # Normalize regression coefficients
    total = sum(reg)
    if total > 0:
        reg = [coef/total for coef in reg]
    
    sample_results = {
        covered_lineages[i]: round(reg[i], 3) for i in range(len(covered_lineages))
    }
    return sample_results, X, Y, covered_muts

