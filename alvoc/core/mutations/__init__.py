from alvoc.core.mutations.prepare import prepare
from alvoc.core.mutations.convert import aa, nt

# Mutation root level functions
def convert_aa(tax_id: str | None, genbank_file: str | None, mut: str, outdir: str) -> list[str]:
    """
    Convert amino acid mutation to nucleotide mutations for a given virus

    Args:
        tax_id : Taxonomic ID of the virus. Required if 'genbank_file' is not provided.
        genbank_file : Path to the GenBank file. Required if 'tax_id' is not provided.
        mut : The amino acid mutation in the format 'GENE:aaPOSITIONaaNEW'.
        outdir : Output directory for results and intermediate data. Defaults to the current directory.

    Returns:
        A list of nucleotide mutations.
    """
    genes, seq = prepare(tax_id, genbank_file, outdir)
    return aa(mut, genes, seq)

def convert_nt(tax_id: str | None, genbank_file: str | None, mut: str, outdir: str) -> str:
    """
    Convert nucleotide mutation to amino acid mutation for a given virus.

    Args:
        tax_id : Taxonomic ID of the virus. Required if 'genbank_file' is not provided.
        genbank_file : Path to the GenBank file. Required if 'tax_id' is not provided.
        mut : The nucleotide mutation in the format 'BASENPOSBASE'.
        outdir : Output directory for results and intermediate data. Defaults to the current directory.

    Returns:
        The amino acid mutation.
    """
    genes, seq = prepare(tax_id, genbank_file, outdir)
    return nt(mut, genes, seq)

def find_mutants(tax_id: str | None, genbank_file: str | None, file_path: str, mutations_path: str | None, min_depth: int, save_img: bool, csv: bool, outdir) -> None:
    """Find mutations in sequencing data, either from BAM files or a sample list. Uses a dictionary of mutation lineages provided as a parameter.

    Args:
        tax_id : Taxonomic ID of the virus. Required if 'genbank_file' is not provided.
        genbank_file : Path to the GenBank file. Required if 'tax_id' is not provided.
        file_path: Path to the file containing sample information or BAM file.
        mutations_path: Path to the file containing mutations or mutation identifier.
        min_depth: Minimum depth for mutation analysis.
        save_img: Whether to save a plot image.
        csv: Whether to generate a CSV file.
        outdir : Output directory for results and intermediate data. Defaults to the current directory.

    Returns:
        Prints number of reads with and without each mutation and generates a heatmap showing their frequencies.
    """
    genes, seq = prepare(tax_id, genbank_file, outdir)
