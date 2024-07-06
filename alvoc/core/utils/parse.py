from alvoc.core.utils.convert import aa


def snv_name(snv: tuple):
    """Generate a name for the SNV based on its components.

    Args:
        snv : A tuple containing the SNV components.

    Returns:
        str: A string representation of the SNV.
    """
    return "{}{}{}".format(*snv)


def parse_mutation(mut: str, genes: dict, seq: str):
    """Parse mutation to handle potential multiple SNVs.

    Args:
        mut : Mutation string, which could include multiple SNVs.
        genes : Dictionary of genes with start and end positions.
        seq : The nucleotide sequence.

    Returns:
        list: A list of parsed SNVs.
    """
    muts = aa(mut, genes, seq) if ":" in mut else [mut]
    return [parse_snv(m) for m in muts]


def parse_mutations(mutations: list[str], genes: dict, seq: str) -> list[str]:
    """
    Parse a list of mutations into single nucleotide changes and amino acid changes.
    Simple nucleotide changes are returned directly, while amino acid changes are
    further processed into nucleotide changes.

    Args:
        mutations : A list of mutation strings.
        genes : Dictionary of genes with start and end positions.
        seq : The nucleotide sequence.
    Returns:
        list: A list of parsed mutations, both nucleotides and processed amino acids.
    """
    # Initialize lists to store nucleotide and amino acid mutations separately
    nts = []
    aas = []

    # Categorize mutations based on their type indicated by the presence of ':'
    for mut in mutations:
        if ":" in mut:
            aas.append(mut)
        else:
            nts.append(mut)

    # Process amino acid mutations to convert them into nucleotide changes
    processed_aas = [
        nt_change for aa_mut in aas for nt_change in aa(aa_mut, genes, seq)
    ]

    # Combine and return the nucleotide changes and processed amino acid changes
    return nts + processed_aas


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
