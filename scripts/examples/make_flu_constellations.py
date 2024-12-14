import aiofiles
import json
from pathlib import Path
import asyncio

async def process_phylogenetic_tree(tree, output_file: Path):
    """
    Process the phylogenetic tree and generate a single JSON file for all clades.
    """
    clades = set()
    traverse_tree(tree["tree"], clades)

    aggregated_data = []

    for clade_name in clades:
        clade_mutations = extract_all_mutations(tree["tree"], clade_name)
        aggregated_data.append(create_clade_entry(clade_name, clade_mutations))

    await write_json_file(aggregated_data, output_file)
    print("Processing completed.")

def traverse_tree(clade, clades):
    """
    Traverse the phylogenetic tree to collect unique clade names.
    """
    if "node_attrs" in clade and "clade_membership" in clade["node_attrs"]:
        clades.add(clade["node_attrs"]["clade_membership"]["value"])

    for child in clade.get("children", []):
        traverse_tree(child, clades)

def extract_all_mutations(clade, target_clade=None):
    """
    Extract all mutations from a specific clade and its parent nodes.
    """
    mutations = set()

    # Collect mutations from this node
    if "branch_attrs" in clade and "mutations" in clade["branch_attrs"]:
        branch_mutations = clade["branch_attrs"]["mutations"]
        if "nuc" in branch_mutations:
            mutations.update(branch_mutations["nuc"])

    # Return if we hit the target clade
    if target_clade and "node_attrs" in clade and "clade_membership" in clade["node_attrs"]:
        clade_membership = clade["node_attrs"]["clade_membership"]["value"]
        if clade_membership == target_clade:
            return list(mutations)

    # Traverse children for additional mutations
    for child in clade.get("children", []):
        child_mutations = extract_all_mutations(child, target_clade)
        mutations.update(child_mutations)

    return list(mutations)

def create_clade_entry(clade_name, clade_mutations):
    """
    Create a data entry for a specific clade.
    """
    entry = {
        "label": clade_name,
        "description": f"{clade_name} defining mutations",
        "sources": [],
        "tags": [clade_name],
        "sites": clade_mutations,
        "note": "Unique mutations for sublineage",
        "rules": {
            "default": {"min_alt": "", "max_ref": ""},
            "Probable": {"min_alt": "", "max_ref": ""},
        },
    }

    return { clade_name : entry}

async def write_json_file(aggregated_data, output_file: Path):
    """
    Write aggregated clade data to a single JSON file.
    """
    async with aiofiles.open(output_file, "w") as file:
        await file.write(json.dumps(aggregated_data, indent=4))
    print(f"Generated JSON file: {output_file}")

async def main(tree_file: str, output_file: Path):
    """
    Main function to process the phylogenetic tree.
    """
    # Ensure the output folder exists
    output_file = Path(output_file)
    output_folder = output_file.parent
    if not output_folder.exists():
        print(f"Creating output directory: {output_folder}")
        output_folder.mkdir(parents=True, exist_ok=True)

    # Read the input tree file
    async with aiofiles.open(tree_file, "r") as file:
        tree_data = json.loads(await file.read())

    # Process the phylogenetic tree and write the JSON file
    await process_phylogenetic_tree(tree_data, output_file)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process Nextstrain phylogenetic tree for clade mutations.")
    parser.add_argument("tree_file", help="Path to the Nextstrain phylogenetic tree JSON file.")
    # parser.add_argument("output_file", help="Path to save the aggregated JSON file.")
    args = parser.parse_args()

    asyncio.run(main(args.tree_file, Path("constellations/flu_constellations.json")))
