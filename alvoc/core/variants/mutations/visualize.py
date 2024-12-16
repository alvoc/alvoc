import matplotlib.pyplot as plt
import seaborn as sns


def plot_mutations(
    df,
    min_depth,
    outdir,
    return_fractions=False,
):
    """Plot mutation fractions across multiple samples using a pandas DataFrame.

    Args:
        df: DataFrame containing mutation results with columns ['mutants', 'sample', 'mutation_count', 'non_mutation_count'].
        min_depth: Minimum read depth to include data in the plot.
        mutants_name: Filename prefix for the saved plot.
        outdir: Base directory for saving the file.
        return_fractions: If True, returns mutation fractions for testing.

    Returns:
        pd.DataFrame: Mutation fractions if return_fractions=True.
    """
    sns.set_theme()

    # Pivot the DataFrame to create a matrix of mutation fractions
    df["total_count"] = df["mutation_count"] + df["non_mutation_count"]
    df["fraction"] = df.apply(
        lambda row: row["mutation_count"] / row["total_count"]
        if row["total_count"] >= min_depth
        else -1,
        axis=1,
    )

    fractions_pivot = df.pivot(index="mutants", columns="sample", values="fraction")

    no_reads = fractions_pivot.map(lambda x: x == -1)

    # Prepare for plotting
    fontsize_pt = plt.rcParams["ytick.labelsize"]
    dpi = 72.27

    matrix_height_pt = fontsize_pt * (fractions_pivot.shape[0] + 30)
    matrix_height_in = matrix_height_pt / dpi

    top_margin = 0.10
    bottom_margin = 0.20
    figure_height = matrix_height_in / (1 - top_margin - bottom_margin)
    figure_width = fractions_pivot.shape[1] * 2 + 5

    fig, ax = plt.subplots(
        figsize=(figure_width, figure_height),
        gridspec_kw=dict(top=1 - top_margin, bottom=bottom_margin),
    )

    sns.heatmap(
        fractions_pivot,
        annot=True,
        mask=no_reads,
        cmap=sns.cm.rocket_r,
        xticklabels=fractions_pivot.columns,
        yticklabels=fractions_pivot.index,
        vmin=0,
        vmax=1,
        fmt=".2f",
    )
    plt.xlabel("Sample")
    plt.xticks(rotation=30, ha="right", rotation_mode="anchor")
    plt.ylabel("Mutation")

    # Save the plot
    img_path = outdir / "mutations.png"
    plt.savefig(img_path, dpi=300)

    if return_fractions:
        return fractions_pivot
