from pathlib import Path

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
    show_stacked: bool = False,
    ts: bool = False,
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
        show_stacked: Whether to show stacked visualizations.
        unique: Whether to consider unique mutations only.
        l2: Whether to use a secondary method for regression analysis.
        ts: Whether to create a time-series plot.

    Returns:
        None: Results are printed or plotted directly.
    """
    # Create or find directory
    out = create_dir(outdir)

    # Extract the genome sequence and gene coordinates for the target virus
    genes, seq = precompute(virus, out)

    # Convert lineage data to mutation-centric format
    mut_lins = parse_lineages(constellations)

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

    sample_results = []
    sample_names = []
    if ".bam" in samples.suffix:
        result = quantify_lineages(
            sample=samples,
            mut_lins=mut_lins,
            genes=genes,
            seq=seq,
            black_list=included_lineages,
            white_list=excluded_lineages,
            min_depth=min_depth,
            unique=unique,
            l2=l2,
            return_data=True,
        )

        if result is None:
            logger.info("No coverage or analysis couldn't be performed.")
        else:
            sr, X, Y, covered_muts = result
            if show_stacked:
                plot_lineage_predictions(sr, X, Y, covered_muts, out)
                plot_lineage_pie(sr, out)
            sample_results.append(sr)
            sample_names.append("")
    else:
        with open(samples, "r") as f:
            sample_df = [line.strip().split("\t") for line in f if line.strip()]
            for bam_path, sample_label in sample_df:
                if bam_path.endswith(".bam"):
                    path = Path(bam_path)
                    sr = quantify_lineages(
                        sample=path,
                        mut_lins=mut_lins,
                        genes=genes,
                        seq=seq,
                        white_list=included_lineages,
                        black_list=excluded_lineages,
                        return_data=False,
                        min_depth=min_depth,
                        unique=unique,
                        l2=l2,
                    )
                    if sr:
                        sample_results.append(sr)
                        sample_names.append(sample_label)

    if sample_results:
        lineage_df = compute_lineage_df(sample_results, sample_names)
        lineage_df.to_csv(out / "lineages.csv", index=False)

        if ts:
            plot_lineages_timeseries(sample_results, sample_names, out)
        else:
            plot_lineages(sample_results, sample_names, out, bool(white_list))


def compute_lineage_df(sample_results, sample_names):
    rows = []
    for i, sample in enumerate(sample_results):
        row = {"Sample name": sample_names[i]}
        for key, value in sample.items():
            row[key] = round(value, 3)
        rows.append(row)
    headers = ["Sample name"] + sorted(
        set(key for sample in sample_results for key in sample.keys())
    )
    df = pd.DataFrame(rows, columns=headers)
    return df


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
    return_data: bool = False,
) -> dict | tuple | None:
    """
    Identify and quantify lineages in a BAM file based on predefined mutations.

    Args:
        bam_path: Path to the BAM file.
        mut_lins: Dictionary containing mutation lineages and their occurrences.
        genes: Dictionary mapping gene names to their start and end positions in the sequence.
        seq: Complete nucleotide sequence as a string.
        black_list: Lineages to ignore.
        return_data: Flag to return detailed data structures.
        min_depth: Minimum read depth to consider a mutation covered.
        lineages: List of lineages to consider.
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
    if return_data:
        return sample_results, X, Y, covered_muts
    return sample_results


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
