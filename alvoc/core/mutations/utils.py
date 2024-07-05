from pysam import PileupColumn
from alvoc.core.mutations.convert import aa

def parse_snv(snv: str):
    """Parse a single nucleotide variant (SNV).

    Args:
        snv : The SNV string in the format 'A123G'.

    Returns:
        tuple: A tuple containing the original base, position, and new base.
    """
    pos = int(snv[1:-1])
    old_bp = snv[0]
    new_bp = snv[-1]
    return old_bp, pos, new_bp


def snv_name(snv: tuple):
    """Generate a name for the SNV based on its components.

    Args:
        snv : A tuple containing the SNV components.

    Returns:
        str: A string representation of the SNV.
    """
    return '{}{}{}'.format(*snv)


def parse_mutation(mut: str, genes: dict, seq: str):
    """Parse mutation to handle potential multiple SNVs.

    Args:
        mut : Mutation string, which could include multiple SNVs.
        genes : Dictionary of genes with start and end positions.
        seq : The nucleotide sequence.

    Returns:
        list: A list of parsed SNVs.
    """
    muts = aa(mut, genes, seq) if ':' in mut else [mut]
    return [parse_snv(m) for m in muts]


def mut_in_col(pileupcolumn: PileupColumn, mut : str):
    """Count the occurrences and non-occurrences of a mutation in a pileup column.

    Args:
        pileupcolumn : A pileup column from a BAM file.
        mut : The mutation to count.

    Returns:
        tuple: A tuple containing counts of mutations and non-mutations.
    """
    muts = not_muts = 0
    for pileupread in pileupcolumn.pileups:
        qpos = pileupread.query_position
        if qpos is None:
            not_muts += 1
            continue
        base = pileupread.alignment.query_sequence
        if base is not None:
            base = base[qpos]
        if base == mut:
            muts += 1
        else:
            not_muts += 1
    return muts, not_muts


def print_mut_results(mut_results, min_depth):
    """Print the results of mutation detection.

    Args:
        mut_results (dict): Dictionary containing mutation results.
        min_depth (int): Minimum depth to consider for reporting.
    """
    cov = mut_cov = 0
    for name, (muts, not_muts) in mut_results.items():
        total = muts + not_muts
        if total >= min_depth:
            cov += 1
            if muts > 0:
                mut_cov += 1

    print('{}/{} mutations covered'.format(cov, len(mut_results)))
    print('{}/{} mutations detected'.format(mut_cov, len(mut_results)))
    
def mut_idx(mut, genes, seq):
    """Get the index of the mutation for sorting purposes.

    Args:
        mut (str): Mutation string.
        genes (dict): Dictionary of genes with start and end positions.
        seq (str): The nucleotide sequence.

    Returns:
        int: The first genomic index of the mutation if available, otherwise -1.
    """
    snvs = parse_mutation(mut, genes, seq)  # Assumes that parse_mutation doesn't require mut_lins
    return snvs[0][1] if snvs else -1

def write_csv(sample_results, sample_names, min_depth, home_lab, mutants_name):
    """Write mutation fractions to a CSV file.

    Args:
        sample_results (list): List of dictionaries containing mutation results.
        sample_names (list): Names of the samples.
        min_depth (int): Minimum read depth to include data in the output.
        home_lab (str): Base directory for saving the file.
        mutants_name (str): Filename suffix for the CSV.
    """
    names = list(sample_results[0].keys())
    csv_headers = ['Mutation'] + [f'{name} %' for name in sample_names]
    csv_rows = [[name] + [str(round(mut_results[name][0] / (mut_results[name][0] + mut_results[name][1]), 4) if (mut_results[name][0] + mut_results[name][1]) >= min_depth else -1) for mut_results in sample_results] for name in names]

    with open(f'{home_lab}_{mutants_name}_mutations.csv', 'w') as file:
        file.write('\n'.join(','.join(row) for row in [csv_headers] + csv_rows))
