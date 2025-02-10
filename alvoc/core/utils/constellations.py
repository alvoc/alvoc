import json
from pathlib import Path
import requests
import os
from collections import defaultdict

def download_phylogenetic_tree(url):
    """
    Download the phylogenetic tree JSON from a specified URL.
    Args:
        url (str): URL to the phylogenetic tree JSON file.
    Returns:
        dict: Parsed JSON data of the phylogenetic tree.
    """
    headers = {
        "Accept": "application/vnd.nextstrain.dataset.main+json"
    }
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    return response.json()

def calculate_clade_mutation_proportions(tree, proportion_threshold=0.90):
    """
    Calculate the proportion of nodes in a clade carrying each mutation.
    Args:
        tree (dict): The phylogenetic tree.
        proportion_threshold (float): Minimum proportion required to include a mutation.
    Returns:
        dict: Dictionary with clade names as keys and sets of significant mutations as values.
    """
    clade_data = defaultdict(set)  # Map clade name -> set of mutations
    clade_node_counts = defaultdict(int)  # Map clade name -> node count
    mutation_counts = defaultdict(lambda: defaultdict(int))  # Map clade name -> mutation -> count

    total_mutations = 0  # Track total mutations across all clades
    retained_mutations = 0  # Track retained mutations after filtering

    def traverse_tree(node, inherited_mutations=set()):
        """
        Recursive function to traverse the tree and propagate mutations.
        Args:
            node (dict): Current node in the tree.
            inherited_mutations (set): Mutations inherited from the parent node.
        """
        # Retrieve the clade name for this node
        clade_name = node.get("node_attrs", {}).get("clade_membership", {}).get("value")
        # Retrieve mutations specific to this branch
        current_mutations = set(node.get("branch_attrs", {}).get("mutations", {}).get("nuc", []))

        # Combine inherited mutations with current mutations
        all_mutations = inherited_mutations | current_mutations

        # Reversions: Handle mutations that restore the root sequence
        # Example: If "A100T" happens in an ancestor and "T100A" happens here, remove it from all_mutations
        for mutation in current_mutations:
            if mutation[1:] + mutation[0] in inherited_mutations:  # Check for a reversion
                all_mutations.remove(mutation)  # Remove reversion

        # Update mutation counts and node counts for this clade
        if clade_name:
            clade_node_counts[clade_name] += 1
            for mutation in all_mutations:
                mutation_counts[clade_name][mutation] += 1

        # Traverse child nodes, passing down the updated mutation set
        for child in node.get("children", []):
            traverse_tree(child, inherited_mutations=all_mutations)

    # Start traversing the tree from the root node
    traverse_tree(tree)

    # Process mutation counts to determine significant mutations for each clade
    original_clades = len(clade_node_counts)
    filtered_clades = 0

    for clade_name, counts in mutation_counts.items():
        total_nodes = clade_node_counts[clade_name]
        for mutation, count in counts.items():
            total_mutations += 1
            if count / total_nodes >= proportion_threshold:
                clade_data[clade_name].add(mutation)
                retained_mutations += 1

        # Track clades with no significant mutations
        if not clade_data[clade_name]:
            filtered_clades += 1

    filtered_mutations = total_mutations - retained_mutations
    remaining_clades = original_clades - filtered_clades

    # Print a summary of the results
    print(
        f"""
        === Summary ===
        Original lineages       : {original_clades}
        Lineages filtered out   : {filtered_clades}
        Lineages remaining      : {remaining_clades}
        ------------------------
        Total mutations         : {total_mutations}
        Mutations filtered out  : {filtered_mutations}
        Mutations retained      : {retained_mutations}
        """
    )

    return clade_data

def create_constellation_entries(clade_data):
    """
    Create constellation entries for each clade in the specified format.
    Args:
        clade_data (dict): Dictionary with clade names and their mutations.
    Returns:
        dict: Dictionary where keys are clade names and values are dictionaries 
              with keys 'lineage', 'label', 'description', 'sources', 'tags', 'sites', and 'note'.
    """
    constellation_entries = {}
    for clade_name, mutations in clade_data.items():
        constellation_entries[clade_name] = {
            "lineage": clade_name,
            "label": f"{clade_name}-like",  # You can adjust this label as needed
            "description": f"{clade_name} lineage defining mutations",
            "sources": [],
            "tags": [clade_name],
            "sites": list(mutations),
            "note": "Unique mutations for sublineage"
        }
    return constellation_entries

def save_constellations_to_json(constellation_entries, output_dir):
    """
    Save constellation entries to a JSON file.
    Args:
        constellation_entries (dict): Dictionary of constellation entries.
        output_dir (str): Directory to save the JSON file.
    """
    output_file = os.path.join(output_dir, "constellations.json")
    with open(output_file, "w") as outfile:
        json.dump(constellation_entries, outfile, indent=4)
    print(f"Constellation JSON file created: {output_file}")

def save_lineages_to_txt(clade_names, output_dir):
    """
    Save clade names to a text file.
    Args:
        clade_names (list): List of clade names.
        output_dir (str): Directory to save the text file.
    """
    output_file = os.path.join(output_dir, "lineages.txt")
    with open(output_file, "w") as outfile:
        for clade_name in sorted(clade_names):
            outfile.write(clade_name + "\n")
    print(f"Lineages file created: {output_file}")

def make_constellations(url : str, outdir: Path, proportion_threshold: int ):
    print("Downloading phylogenetic tree...")
    tree_data = download_phylogenetic_tree(url)

    print("Calculating clade mutation proportions...")
    tree = tree_data["tree"]  # Assuming the tree starts here
    clade_data = calculate_clade_mutation_proportions(tree, proportion_threshold=proportion_threshold)

    # Remove clades with no significant mutations
    empty_lineages = [lineage for lineage, mutations in clade_data.items() if not mutations]
    for lineage in empty_lineages:
        del clade_data[lineage]

    print("Creating constellation entries...")
    constellation_entries = create_constellation_entries(clade_data)

    print("Saving constellation JSON...")
    save_constellations_to_json(constellation_entries, outdir)

    print("Saving lineages to text file...")
    save_lineages_to_txt(clade_data.keys(), outdir)

    print("Processing complete.")