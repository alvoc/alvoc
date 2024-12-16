from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from ortools.linear_solver import pywraplp
from sklearn.linear_model import LinearRegression

from alvoc.core.utils import create_dir, logging
from alvoc.core.lineages.prepare import parse_lineages
from alvoc.core.lineages.visualize import (
    plot_lineages,
    plot_lineages_timeseries,
    plot_lineage_pie,
    plot_lineage_predictions,
)
from alvoc.core.mutations.analyze import find_mutants_in_bam
from alvoc.core.utils.parse import parse_mutations
from alvoc.core.utils.precompute import precompute

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
        process_sample(
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

    # Save results to a CSV file
    results_df.to_csv(out / "lineages.csv", index=False)
    # plot_lineages(sample_results, sample_names, out, bool(white_list))


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

    sample_results = {
        covered_lineages[i]: round(reg[i], 3) for i in range(len(covered_lineages))
    }
    return sample_results, X, Y, covered_muts


def do_regression(lmps: list[list[float]], Y: np.ndarray):
    """
    Perform linear regression on lineage mutation profiles against observed data,
    adjusting the model iteratively if the sum of coefficients exceeds 1 to find the best model
    based on regression score that satisfies this constraint.

    Args:
        lmps: List of lineage mutation profiles.
        Y: Observed mutation frequencies.

    Returns:
        tuple: Tuple containing the matrix of mutation profiles and the best regression coefficients.
    """
    X = np.array(lmps).T
    reg = LinearRegression(fit_intercept=False, positive=True).fit(X, Y)
    best_score = reg.score(X, Y)
    best_reg = reg

    if sum(reg.coef_) > 1:
        for i in range(len(X)):
            if Y[i] == 0:
                new_X = np.concatenate((X[:i], X[i + 1 :]))
                new_Y = np.concatenate((Y[:i], Y[i + 1 :]))
                new_reg = LinearRegression(fit_intercept=False, positive=True).fit(
                    new_X, new_Y
                )
                new_score = new_reg.score(new_X, new_Y)
                if sum(new_reg.coef_) <= 1 and new_score > best_score:
                    best_reg = new_reg
                    best_score = new_score

    return X, [coef for coef in best_reg.coef_]


def do_regression_linear(lmps, Y, muts):
    """
    Perform a linear programming approach to minimize the error between predicted and observed mutation frequencies.

    Args:
        lmps (list of list of floats): Lineage mutation profiles.
        Y (list of floats): Observed mutation frequencies.
        muts (list of str): List of mutation identifiers.

    Returns:
        tuple: Matrix of mutation profiles, solution values for lineages, and a dictionary of mutation differences.
    """

    solver = pywraplp.Solver.CreateSolver("GLOP")
    num_lins = len(lmps)
    num_muts = len(lmps[0])
    lins = [solver.NumVar(0, 1, f"x_{i}") for i in range(num_lins)]
    t = [solver.NumVar(0, solver.infinity(), f"t_{i}") for i in range(num_muts)]

    # Set up constraints
    for i, mut_freq in enumerate(Y):
        constraint_less = solver.Constraint(-solver.infinity(), mut_freq)
        constraint_more = solver.Constraint(mut_freq, solver.infinity())
        for j in range(num_lins):
            constraint_less.SetCoefficient(lins[j], lmps[j][i])
            constraint_more.SetCoefficient(lins[j], lmps[j][i])
        constraint_less.SetCoefficient(t[i], -1)
        constraint_more.SetCoefficient(t[i], 1)

    # Objective: minimize the sum of t, which represents the residuals
    objective = solver.Objective()
    for ti in t:
        objective.SetCoefficient(ti, 1)
    objective.SetMinimization()

    solver.Solve()

    X = np.array(lmps).T
    mut_diffs = {
        mut: Y[i] - sum(lins[j].solution_value() * lmps[j][i] for j in range(num_lins))
        for i, mut in enumerate(muts)
    }

    return X, [lin.solution_value() for lin in lins], mut_diffs
