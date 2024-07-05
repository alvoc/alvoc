from alvoc.core.utils.convert import aa

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

def mut_idx(mut, genes, seq) -> int:
    """Get the index of the mutation for sorting purposes.

    Args:
        mut (str): Mutation string.
        genes (dict): Dictionary of genes with start and end positions.
        seq (str): The nucleotide sequence.

    Returns:
        int: The first genomic index of the mutation if available, otherwise -1.
    """
    snvs = parse_mutation(mut, genes, seq)
    return snvs[0][1] if snvs else -1