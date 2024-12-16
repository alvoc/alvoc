from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

sns.set_theme()


def plot_depths(results_df: pd.DataFrame, inserts: list[list], outdir: Path):
    """
    Plots the logarithmic depth of sequencing results across samples and pools using a bar plot.

    Args:
        results_df: Pandas DataFrame containing `sample`, `amplicon_id`, and `depth` columns.
        inserts: A list of data used for generating additional metadata about the amplicons.
        outdir: Output directory for plot.
    """
    # Add "Pool" information based on amplicon number (example logic: even = Pool 1, odd = Pool 2)
    results_df["Pool"] = results_df["amplicon_id"].apply(
        lambda x: "Pool 1" if int(x.split("_")[-1]) % 2 == 0 else "Pool 2"
    )
    # Add log-transformed depth
    results_df["Log depth"] = results_df["depth"].apply(lambda x: np.log(max(x, 1)))

    # Create a bar plot of log depth by amplicon number
    grid = sns.FacetGrid(
        results_df,
        row="sample",
        hue="Pool",
        height=1.7,
        aspect=8,
        sharey=False,
    )
    grid.map(
        sns.barplot,
        "amplicon_id",
        "Log depth",
        order=sorted(
            results_df["amplicon_id"].unique(), key=lambda x: int(x.split("_")[-1])
        ),
        hue_order=["Pool 1", "Pool 2"],
    )
    grid.add_legend()
    plt.xticks(rotation=90)
    plt.savefig(outdir / "depths.png", dpi=300)


def plot_depths_gc(results_df: pd.DataFrame, outdir: Path):
    """
    Plots the relationship between GC content and log depth for each sample using regression plots,
    and computes the Pearson correlation coefficient.

    Args:
        results_df: Pandas DataFrame containing `sample`, `amplicon_id`, `depth`, and `gc_content` columns.
        outdir: Output directory for plot.
    """
    # Add log-transformed depth to DataFrame
    results_df["Log depth"] = results_df["depth"].apply(lambda x: np.log(max(x, 1)))

    # Ensure there is enough data for plotting
    if results_df.empty:
        raise ValueError(
            "No data available for plotting. Check the `results_df` DataFrame."
        )

    # Create regression plots for GC content vs log depth
    grid = sns.FacetGrid(results_df, col="sample", height=4, aspect=1)
    grid.map(sns.regplot, "gc_content", "Log depth")

    # Compute Pearson correlation for GC content vs log depth
    for sample, group in results_df.groupby("sample"):
        gcs = group["gc_content"]
        depths = group["Log depth"]
        if len(gcs) > 1:
            statistic, pvalue = stats.pearsonr(gcs, depths)
            print(f"Sample: {sample}")
            print(f"  Pearson correlation statistic: {statistic:.3f}")
            print(f"  Pearson p-value: {pvalue:.3g}")
        else:
            print(
                f"Sample: {sample} - Insufficient data to compute Pearson correlation."
            )

    # Save the plot
    plt.savefig(outdir / "gc_depths.png", dpi=300)
