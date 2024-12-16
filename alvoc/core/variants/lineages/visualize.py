from datetime import date
from math import ceil
from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from numpy import ndarray
import logging

logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)


def plot_lineages(
    data: pd.DataFrame,
    outdir: Union[str, Path],
    all_lins: bool = False,
    xlabel: str = "Samples",
    ylabel: str = "Lineage",
) -> None:
    """
    Plot a heatmap of lineage distributions per sample.

    Args:
        data: DataFrame with columns ['sample', 'lineage', 'abundance'].
        outdir: Output directory for the plot.
        all_lins: Flag to show all lineages or only significant ones.
    """
    sns.set_theme()

    # Filter lineages based on the `all_lins` flag
    if not all_lins:
        data = data[data["abundance"] > 0.01]

    # Pivot the data for heatmap compatibility
    pivot_table = (
        data.pivot(index="lineage", columns="sample", values="abundance").fillna(0)
        * 100
    )
    sample_names = pivot_table.columns.tolist()
    lineage_names = pivot_table.index.tolist()
    lin_fractions = pivot_table.values

    # Calculate figure dimensions dynamically
    fontsize_pt = plt.rcParams["ytick.labelsize"]
    dpi = 72.27
    matrix_height_pt = fontsize_pt * (len(lineage_names) + 30)
    matrix_height_in = matrix_height_pt / dpi
    top_margin = 0.10
    bottom_margin = 0.20
    figure_height = matrix_height_in / (1 - top_margin - bottom_margin)
    figure_width = len(sample_names) * 2 + 5

    # Create the heatmap
    fig, ax = plt.subplots(
        figsize=(figure_width, figure_height),
        gridspec_kw={"top": 1 - top_margin, "bottom": bottom_margin},
    )
    _ = sns.heatmap(
        lin_fractions,
        annot=True,
        cmap=sns.cm.rocket_r,
        xticklabels=sample_names,
        yticklabels=lineage_names,
        vmin=0,
        vmax=100,
        cbar_kws={"format": "%.0f%%"},
        fmt=".1f",
    )

    # Configure labels and ticks
    plt.xlabel(xlabel)
    plt.xticks(rotation=30, ha="right", rotation_mode="anchor")
    plt.ylabel(ylabel)
    plt.subplots_adjust(bottom=0.3, left=0.6)

    # Save the plot
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_file = outdir / "lineage_abundance_heatmap.png"
    plt.savefig(output_file, dpi=300)
    plt.close()


def plot_lineages_timeseries(
    sample_results: list[dict], sample_names: list[str], outdir: Path
):
    """Plot time series data for lineages across multiple samples.

    Args:
        sample_results: List of dictionaries with sample results.
        sample_names: List of sample names including date and location tags.
        outdir : Output directory for plot.

    """
    sns.set_theme()

    # Determine unique locations and corresponding dates
    locations, dates = [], []
    for sample_name in sample_names:
        loc, dt = sample_name.split("_")
        dates.append(date.fromisoformat(dt))
        locations.append(loc)

    loc_set = sorted(set(locations))
    ncols = 3
    nrows = ceil(len(loc_set) / ncols)
    fig, axes = plt.subplots(nrows, ncols)
    plt.tight_layout()

    # Prepare and plot data for each location
    for i, loc in enumerate(loc_set):
        r, c = divmod(i, ncols)
        loc_results = [sr for j, sr in enumerate(sample_results) if locations[j] == loc]
        loc_dates = [dt for j, dt in enumerate(dates) if locations[j] == loc]
        data = {
            name: [sr.get(name, 0) for sr in loc_results]
            for name in sorted(sample_results[0].keys())
        }
        df = pd.DataFrame(data, index=loc_dates)
        axes[r, c].set_title(loc)
        df.plot.area(ax=axes[r, c], fontsize=8, legend=False, rot=30)
        if c == ncols - 1 or i == len(loc_set) - 1:
            axes[r, c].legend(loc="upper left")

    plt.savefig(outdir / "lineages_timeseries.png", dpi=300)


def plot_lineage_predictions(
    sample_results: dict, X: ndarray, Y: ndarray, covered_muts: list, outdir: Path
):
    """Display predicted versus observed lineage contributions for mutations.

    Args:
        sample_results: Lineage contribution results.
        X: Lineage mutation profiles.
        Y: Observed mutation frequencies.
        covered_muts (list): Mutations that meet coverage requirements.
        outdir : Output directory for plot.
    """
    sns.set_theme()

    merged_lins = list(sample_results.keys())
    d = {ml: [] for ml in merged_lins}
    d["Observed"] = []
    idx = []

    # Compute predictions and compare with observed data
    for i, cm in enumerate(covered_muts):
        Y_p = sum(X[i][j] * sample_results[ml] for j, ml in enumerate(merged_lins))
        if not (Y[i] < 0.02 and Y_p < 0.02):
            idx.extend([f"{cm}\nobserved", f"{cm}\npredicted"])
            d["Observed"].extend([Y[i], 0])
            for j, ml in enumerate(merged_lins):
                d[ml].extend([0, X[i][j] * sample_results[ml]])

    df = pd.DataFrame(d, index=idx)
    df.plot.bar(stacked=True, rot=0)
    plt.tight_layout()
    plt.savefig(outdir / "lineage_predictions.png", dpi=300)


def plot_lineage_pie(sample_results: dict, outdir: Path):
    """Plot a pie chart of lineage distributions.

    Args:
        sample_results: Lineage distribution data.
        outdir : Output directory for plot.
    """
    sns.set_theme()

    merged_lins = list(sample_results.keys())
    freqs = [sample_results[ml] for ml in merged_lins]
    freqs = [f / sum(freqs) for f in freqs]  # Normalize frequencies
    df = pd.DataFrame({"Fraction": freqs}, index=merged_lins)
    df = df[df["Fraction"] > 0]  # Filter zero fractions
    df.plot.pie(y="Fraction", legend=False, autopct="%1.1f%%", ylabel="")
    plt.savefig(outdir / "lineage_pie.png", dpi=300)
