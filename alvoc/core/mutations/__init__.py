from alvoc.core.mutations.prepare import prepare
from alvoc.core.mutations.convert import aa, nt
from alvoc.core.mutations.analyze import find_mutants as find_mutants_internal

# Mutation root level functions
def convert_aa(tax_id: str | None, genbank_file: str | None, mut: str, outdir: str):
    """
    Convert amino acid mutation to nucleotide mutations for a given virus

    Args:
        tax_id (str, optional): Taxonomic ID of the virus. Required if 'genbank_file' is not provided.
        genbank_file (str, optional): Path to the GenBank file. Required if 'tax_id' is not provided.
        mut (str): The amino acid mutation in the format 'GENE:aaPOSITIONaaNEW'.
        outdir (str, optional): Output directory for results and intermediate data. Defaults to the current directory.

    Returns:
        list: A list of nucleotide mutations.
    """
    genes, seq = prepare(tax_id, genbank_file, outdir)
    return aa(mut, genes, seq)

def convert_nt(tax_id: str | None, genbank_file: str | None, mut: str, outdir: str):
    """
    Convert nucleotide mutation to amino acid mutation for a given virus.

    Args:
        tax_id (str, optional): Taxonomic ID of the virus. Required if 'genbank_file' is not provided.
        genbank_file (str, optional): Path to the GenBank file. Required if 'tax_id' is not provided.
        mut (str): The nucleotide mutation in the format 'BASENPOSBASE'.
        outdir (str, optional): Output directory for results and intermediate data. Defaults to the current directory.

    Returns:
        str: The amino acid mutation.
    """
    genes, seq = prepare(tax_id, genbank_file, outdir)
    return nt(mut, genes, seq)
