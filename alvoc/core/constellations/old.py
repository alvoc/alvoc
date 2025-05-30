import json
from pathlib import Path
import requests
import os
from collections import defaultdict, deque

def download_phylogenetic_tree(url):
    """
    Download the phylogenetic tree JSON from a specified URL.
    Args:
        url (str): URL to the phylogenetic tree JSON file.
    Returns:
        dict: Parsed JSON data of the phylogenetic tree.
    """
    headers = {"Accept": "application/vnd.nextstrain.dataset.main+json"}
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    return response.json()


def extract_raw_mutation_counts(tree):
    """
    Walk the Nextstrain JSON tree *iteratively* and count every NA mutation per clade
    using a bitmask for each mutation‐set. Orders of magnitude faster than Python set unions.

    Returns:
      - clade_node_counts: dict str→int
      - mutation_counts : dict str→list[int]  (same order as `mutations_list`)
      - mutations_list  : list[str]            (all unique mutations encountered)
    """

    # 1) First pass: collect every mutation string into a list
    mutations = set()
    q = deque([tree])
    while q:
        node = q.popleft()
        muts = node.get("branch_attrs", {}).get("mutations", {}).get("nuc", [])
        mutations.update(muts)
        q.extend(node.get("children", []))
    mutations_list = sorted(mutations)
    idx_map = {m: i for i, m in enumerate(mutations_list)}

    # 2) Prepare counters
    clade_node_counts = defaultdict(int)
    # per‐clade list of counts per mutation index
    mutation_counts = defaultdict(lambda: [0] * len(mutations_list))

    # 3) Traverse tree iteratively, carrying a bitmask of inherited mutations
    q = deque([(tree, 0)])  # (node, inherited_mask)
    while q:
        node, inh_mask = q.popleft()

        # build mask for this branch’s nukes
        curr_mask = 0
        for m in node.get("branch_attrs", {}).get("mutations", {}).get("nuc", []):
            i = idx_map[m]
            curr_mask |= 1 << i

        # union inherited + current
        all_mask = inh_mask | curr_mask

        # === drop true reversions ===
        tmp = curr_mask
        while tmp:
            lowbit = tmp & -tmp
            i = lowbit.bit_length() - 1
            mut = mutations_list[i]
            # parse "A100T" → ref='A', pos='100', alt='T'
            ref, pos, alt = mut[0], mut[1:-1], mut[-1]
            reverse = f"{alt}{pos}{ref}"
            if reverse in idx_map:
                j = idx_map[reverse]
                if (inh_mask >> j) & 1:            # if reverse was inherited
                    all_mask &= ~(1 << i)         # drop this forward mutation
            tmp ^= lowbit
        # ================================


        # update counts for this clade
        clade = node.get("node_attrs", {}).get("clade_membership", {}).get("value")
        if clade:
            clade_node_counts[clade] += 1
            # iterate bits set in all_mask
            mask = all_mask
            while mask:
                lowbit = mask & -mask
                bit_i = lowbit.bit_length() - 1
                mutation_counts[clade][bit_i] += 1
                mask ^= lowbit

        # enqueue children
        for child in node.get("children", []):
            q.append((child, all_mask))

    return clade_node_counts, mutation_counts, mutations_list


def filter_by_proportion(clade_node_counts, mutation_counts, mutations_list, threshold):
    """
    Replace your old filter_by_proportion.
    Prints the same summary, then returns profiles: dict clade→set(mutation_str).
    """
    original_clades = len(clade_node_counts)
    total_mutations = sum(len(cnts) for cnts in mutation_counts.values())

    profiles = {}
    for clade, counts in mutation_counts.items():
        n = clade_node_counts[clade]
        # pick mutation indices where freq >= threshold
        kept = {mutations_list[i] for i, c in enumerate(counts) if c / n >= threshold}
        if kept:
            profiles[clade] = kept

    filtered_clades = original_clades - len(profiles)
    retained_mutations = sum(len(s) for s in profiles.values())
    filtered_mutations = total_mutations - retained_mutations
    remaining_clades = len(profiles)

    print(f"""
    === Summary ===
    Original lineages       : {original_clades}
    Lineages filtered out   : {filtered_clades}
    Lineages remaining      : {remaining_clades}
    ------------------------
    Total mutations         : {total_mutations}
    Mutations filtered out  : {filtered_mutations}
    Mutations retained      : {retained_mutations}
    """)

    return profiles


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
            "note": "Unique mutations for sublineage",
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


def make_constellations(
    url: str,
    outdir: Path,
    proportion_threshold: float,
):
    # 1) Fetch the Nextstrain JSON
    print("Downloading phylogenetic tree…")
    tree_data = download_phylogenetic_tree(url)
    tree = tree_data["tree"]

    # 2) Walk the tree and count every mutation per clade
    print("Extracting raw mutation counts…")
    node_counts, mutation_counts, mutations_list = extract_raw_mutation_counts(tree)

    # 3) Filter by your proportion_threshold (prints its own summary)
    print(
        f"Filtering mutations to those fixed in ≥{proportion_threshold*100:.0f}% of nodes…"
    )
    profiles = filter_by_proportion(
        clade_node_counts=node_counts,
        mutation_counts=mutation_counts,
        mutations_list=mutations_list,
        threshold=proportion_threshold,
    )

    # 4) Build the constellation entries
    constellation_entries = create_constellation_entries(profiles)
     
    # 5) Write out the files
    print("Saving constellation JSON…")
    save_constellations_to_json(constellation_entries, outdir)
    print("Saving lineage list…")
    save_lineages_to_txt(constellation_entries.keys(), outdir)

    print("Processing complete.")
