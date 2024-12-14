def parse_lineages(lineage_data : dict):
    """
    Convert lineage-centric JSON format to mutation-centric format with 0s for absent mutations.

    Args:
        lineage_data (dict): Original lineage data format
        output_file (Path): Path to save the mutation-centric JSON file.
    """
    # Gather all unique mutations and lineages
    all_mutations = set()
    all_lineages = lineage_data.keys()

    for lineage, details in lineage_data.items():
        all_mutations.update(details["sites"])

    # Create mutation-centric dictionary
    mutation_centric_data = {mutation: {} for mutation in all_mutations}

    for mutation in all_mutations:
        for lineage in all_lineages:
            lineage_label = f"{lineage}-like"
            if mutation in lineage_data[lineage]["sites"]:
                mutation_centric_data[mutation][lineage_label] = 1
            else:
                mutation_centric_data[mutation][lineage_label] = 0

    return mutation_centric_data