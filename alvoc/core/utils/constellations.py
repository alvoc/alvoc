import json
from pathlib import Path
import requests
import os
from collections import Counter, defaultdict, deque
import sys
import logging


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


def prune_lineages_without_unique_sites(constellation_entries, idx_map, min_unique: int = 1):
    """
    Return a new dict of constellation_entries containing only those lineages
    that have ≥ min_unique private sites.  The others are silently dropped.
    """
    logger = logging.getLogger(__name__)

    # 1) Build per-lineage bitmask
    clade_masks = {}
    for lin, cfg in constellation_entries.items():
        m = 0
        for mut in cfg["sites"]:
            m |= 1 << idx_map[mut]
        clade_masks[lin] = m

    # 2) Count in how many lineages each bit (site) appears
    idx_counts = Counter()
    for mask in clade_masks.values():
        tmp = mask
        while tmp:
            lowbit = tmp & -tmp
            bit_i = lowbit.bit_length() - 1
            idx_counts[bit_i] += 1
            tmp ^= lowbit

    # 3) Build a mask of “unique” sites (count == 1)
    unique_mask = 0
    for bit_i, cnt in idx_counts.items():
        if cnt == 1:
            unique_mask |= 1 << bit_i

    # 4) Keep only those lineages that pass
    pruned = {}
    for lin, mask in clade_masks.items():
        if (mask & unique_mask).bit_count() >= min_unique:
            pruned[lin] = constellation_entries[lin]
        else:
            logger.warning("Dropping lineage %s: only %d unique sites",
                           lin, (mask & unique_mask).bit_count())

    logger.info("Pruned %d → %d lineages (min_unique=%d)",
                len(constellation_entries), len(pruned), min_unique)
    return pruned

def require_unique_sites(constellation_entries, min_unique: int = 1):
    """
    Ensure each lineage has at least `min_unique` unique sites.
    Exits with error if any lineage has fewer.
    """
    logger = logging.getLogger(__name__)

    # 1) Build a map lineage -> list of sites
    profiles = { lin: cfg["sites"] 
                 for lin, cfg in constellation_entries.items() }

    # 2) Count in how many lineages each site appears
    site_counts = Counter()
    for sites in profiles.values():
        # each lineage only contributes once per site
        for s in set(sites):
            site_counts[s] += 1

    # 3) For each lineage, check if it has ≥ min_unique private sites,
    #    bailing out as soon as we hit our quota.
    bad = []
    for lin, sites in profiles.items():
        if min_unique == 1:
            # fast any() pattern
            has_private = any(site_counts[s] == 1 for s in sites)
            if not has_private:
                bad.append((lin, 0))
        else:
            # count up to min_unique
            cnt = 0
            for s in sites:
                if site_counts[s] == 1:
                    cnt += 1
                    if cnt >= min_unique:
                        break
            if cnt < min_unique:
                bad.append((lin, cnt))

    if bad:
        msg = "\n".join(f"  • {lin}: {n} unique sites" for lin, n in bad)
        sys.stderr.write(
            f"ERROR: the following lineages have fewer than {min_unique} unique sites:\n"
            f"{msg}\n\n"
            "Adjust --min-unique-sites or enrich your constellation.\n"
        )
        sys.exit(1)

    logger.info("All %d lineages have ≥%d unique sites",
                len(profiles), min_unique)

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
    min_unique: int,
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
    
    # # 5) Now check those entries
    # print(f"Ensuring each lineage has at least {min_unique} unique site(s)…")
    # idx_map = { m: i for i, m in enumerate(mutations_list) }
    # constellation_entries = prune_lineages_without_unique_sites(constellation_entries, idx_map, min_unique)
    
    # 6) Write out the files
    print("Saving constellation JSON…")
    save_constellations_to_json(constellation_entries, outdir)
    print("Saving lineage list…")
    save_lineages_to_txt(constellation_entries.keys(), outdir)

    print("Processing complete.")
