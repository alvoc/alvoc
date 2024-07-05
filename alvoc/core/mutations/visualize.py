import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def plot_mutations(sample_results, sample_names, min_depth, mutants_name, outdir, return_fractions = False):
    """Plot mutation fractions across multiple samples.

    Args:
        sample_results (list): List of dictionaries containing mutation results.
        sample_names (list): Names of the samples.
        min_depth (int): Minimum read depth to include data in the plot.
        mutants_name (str): Filename prefix for the CSV.
        outdir (str): Base directory for saving the file.
        return_fractions (bool): test parameter for retrieving fractions
    """
    sns.set_theme()
    names = list(sample_results[0].keys())
    sample_counts = [
        [mut_results[mut] for mut in names] for mut_results in sample_results
    ]
    num_mutations = len(names) 

    mut_fractions = [[] for _ in range(num_mutations)]
    for i in range(num_mutations):
        for counts in sample_counts:
            count = counts[i]
            total = count[0] + count[1]
            fraction = count[0]/total if total >= min_depth else -1
            mut_fractions[i].append(round(fraction, 4))

    if return_fractions:
        return mut_fractions  # Return the fractions for testing

    no_reads = np.array([[f == -1 for f in fractions] for fractions in mut_fractions])

    # get the tick label font size
    fontsize_pt = plt.rcParams["ytick.labelsize"]
    dpi = 72.27

    # comput the matrix height in points and inches
    matrix_height_pt = fontsize_pt * (len(mut_fractions) + 30)
    matrix_height_in = matrix_height_pt / dpi

    # compute the required figure height
    top_margin = 0.10  # in percentage of the figure height
    bottom_margin = 0.20  # in percentage of the figure height
    figure_height = matrix_height_in / (1 - top_margin - bottom_margin)
    figure_width = len(sample_names) * 2 + 5

    # build the figure instance with the desired height
    fig, ax = plt.subplots(
        figsize=(figure_width, figure_height),
        gridspec_kw=dict(top=1 - top_margin, bottom=bottom_margin),
    )

    sns.heatmap(
        mut_fractions,
        annot=True, 
        mask=no_reads,
        cmap=sns.cm.rocket_r,
        xticklabels=sample_names,
        yticklabels=names,
        vmin=0,
        vmax=1,
        fmt=".2f",
    )
    plt.xlabel("Sample")
    plt.xticks(rotation=30, ha="right", rotation_mode="anchor")
    plt.ylabel("Mutation") 
    img_path = outdir / f"{mutants_name}_mutations.png"
    if img_path is not None:
        plt.savefig(img_path, dpi=300)
    else:
        plt.show()
   