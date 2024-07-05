def write_csv(sample_results, sample_names, min_depth, mutants_name, outdir):
    """Write mutation fractions to a CSV file.

    Args:
        sample_results (list): List of dictionaries containing mutation results.
        sample_names (list): Names of the samples.
        min_depth (int): Minimum read depth to include data in the output.
        mutants_name (str): Filename prefix for the CSV.
        outdir (str): Base directory for saving the file.
    """
    names = list(sample_results[0].keys())
    csv_headers = ['Mutation'] + [f'{name} %' for name in sample_names]
    csv_rows = [[name] + [str(round(mut_results[name][0] / (mut_results[name][0] + mut_results[name][1]), 4) if (mut_results[name][0] + mut_results[name][1]) >= min_depth else -1) for mut_results in sample_results] for name in names]
    csv_path = outdir / f"{mutants_name}_mutations.csv"

    with open(csv_path, 'w') as file:
        file.write('\n'.join(','.join(row) for row in [csv_headers] + csv_rows))
