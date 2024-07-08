from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

sns.set_theme()


def plot_depths(
    sample_results: list[dict[str, str]],
    sample_names: list[str],
    inserts: list,
    outdir: Path,
):
    """
    Plots the logarithmic depth of sequencing results across samples and pools using a bar plot.

    Args:
    sample_results (list): A list of dictionaries where keys are amplicon numbers and values are depths.
    sample_names (list): A list of names corresponding to the samples in `sample_results`.
    inserts (list): A list of data used for generating additional sample metadata.
    outdir : Output directory for plot.

    """
    samples = sum([len(inserts) * [name] for name in sample_names], [])
    amplicons = sum([list(amplicon.keys()) for amplicon in sample_results], [])
    pools = ["Pool 1" if True else "Pool 2" for _ in amplicons]
    depths = sum(
        [[np.log(max(a, 1)) for a in amplicon.values()] for amplicon in sample_results],
        [],
    )
    data = {
        "Sample": samples,
        "Amplicon number": amplicons,
        "Pool": pools,
        "Log depth": depths,
    }
    df = pd.DataFrame(data)
    grid = sns.FacetGrid(df, row="Sample", hue="Pool", height=1.7, aspect=8)
    grid.map(
        sns.barplot,
        "Amplicon number",
        "Log depth",
        order=[str(i) for i in range(1, 99)],
        hue_order=["Pool 1", "Pool 2"],
    )
    plt.locator_params(axis="x")
    plt.savefig(outdir / "depths.png", dpi=300)
    plt.show()


def plot_depths_gc(
    sample_results: list[dict[str, str]],
    sample_names: list[str],
    inserts: list,
    outdir: Path
):
    """
    Plots the relationship between GC content and log depth for each sample using regression plots,
    and computes the Pearson correlation coefficient.

    Args:
    sample_results (list): A list of dictionaries where keys are amplicon numbers and values are depths.
    sample_names (list): A list of names corresponding to the samples in `sample_results`.
    inserts (list): A list of data used for extracting GC content.
    outdir : Output directory for plot.

    """
    depths = sample_results[0]
    samples = sum([98 * [name] for name in sample_names], [])
    amplicons = sum([list(amplicon.keys()) for amplicon in sample_results], [])
    gcs = [inserts[int(a) - 1][6] for a in amplicons]
    pools = ["Pool 1" if int(a) % 2 == 0 else "Pool 2" for a in amplicons]
    depths = sum(
        [[np.log(max(a, 1)) for a in amplicon.values()] for amplicon in sample_results],
        [],
    )
    data = {
        "Sample": samples,
        "Amplicon number": [int(a) for a in amplicons],
        "GC content": gcs,
        "Pool": pools,
        "Log depth": depths,
    }
    df = pd.DataFrame(data)
    grid = sns.FacetGrid(df, col="Sample")
    grid.map(sns.regplot, "GC content", "Log depth")
    statistic, pvalue = stats.pearsonr(gcs, depths)
    print(f"Pearson correlation statistic: {statistic}")
    print(f"Pearson p-value (probability that variables are independent): {pvalue}")
    plt.savefig(outdir / "gc_depths.png", dpi=300)
    plt.show()
