import typer

from alvoc.core.logging import init_logger
from alvoc.core.mutations.analyze import find_mutants as find_mutants_internal
from alvoc.core.precompute import precompute
from alvoc.core.utils.convert import aa, nt

# from alvoc.core.lineages import find_lineages
# from alvoc.core.amplicons import amplicon_coverage, gc_depth

cli = typer.Typer(
    help="Identify frequencies of concerning mutations from aligned reads"
)


@cli.callback()
def callback():
    init_logger()

    pass


@cli.command()
def convert_amino_acid(
    tax_id=typer.Option(
        None, help="Taxonomic ID of the virus, required if no GenBank file is provided"
    ),
    genbank_file=typer.Option(
        None, help="Path to the GenBank file, required if no tax ID is provided"
    ),
    mut: str = typer.Option(
        ..., help="Amino acid mutation in the format 'GENE:aaPOSITIONaaNEW'"
    ),
    outdir: str = typer.Option(
        ".", help="Output directory for results and intermediate data"
    ),
):
    """
    Convert amino acid mutation to nucleotide mutations for a given virus.
    """
    genes, seq, _ = precompute(tax_id, genbank_file, outdir)
    return aa(mut, genes, seq)


@cli.command()
def convert_nucleotide(
    tax_id=typer.Option(
        None, help="Taxonomic ID of the virus, required if no GenBank file is provided"
    ),
    genbank_file=typer.Option(
        None, help="Path to the GenBank file, required if no tax ID is provided"
    ),
    mut: str = typer.Option(
        ..., help="Nucleotide mutation in the format 'BASENPOSBASE'"
    ),
    outdir: str = typer.Option(
        ".", help="Output directory for results and intermediate data"
    ),
):
    """
    Convert nucleotide mutation to amino acid mutation for a given virus.
    """
    genes, seq, _ = precompute(tax_id, genbank_file, outdir)
    return nt(mut, genes, seq)


@cli.command()
def find_mutants(
    samples_path: str,
    tax_id: str = typer.Option(
        None, help="Taxonomic ID of the virus, required if no GenBank file is provided"
    ),
    genbank_file: str = typer.Option(
        None, help="Path to the GenBank file, required if no tax ID is provided"
    ),
    mutations_path: str = typer.Option(
        None, "--mutations-path", "-m", help="Path to mutations"
    ),
    min_depth: int = typer.Option(40, "--min-depth", "-d", help="Minimum depth"),
    outdir: str = typer.Option(
        ".", help="Output directory for results and intermediate data"
    ),
):
    """
    Find mutations in sequencing data, either from BAM files or a sample list. Uses a dictionary of mutation lineages provided as a parameter.
    """
    genes, seq, out = precompute(tax_id, genbank_file, outdir)
    mut_lins = {}
    find_mutants_internal(
        samples_path, mutations_path, min_depth, mut_lins, genes, seq, out
    )


# @cli.command()
# def find_lineages_command(
#     samples_path: str,
#     lineages_path: str = typer.Option(None, "--lineages-path", "-l", help='Path to lineages'),
#     ts: bool = typer.Option(False, "--ts", "-t", help='Time series'),
#     csv: bool = typer.Option(False, "--csv", "-c", help='Save as CSV'),
#     min_depth: int = typer.Option(40, "--min-depth", "-d", help='Minimum depth'),
#     show_stacked: bool = typer.Option(False, "--show-stacked", "-s", help='Show stacked'),
#     unique: bool = typer.Option(False, "--unique", "-u", help='Unique'),
#     save_img: bool = typer.Option(False, "--save-img", "-s", help='Save image'),
#     l2: bool = typer.Option(False, "--l2", "-l2", help='Level 2')
# ):
#     """Find lineages in samples"""
#     find_lineages(samples_path, lineages_path=lineages_path, ts=ts, csv=csv, min_depth=min_depth, show_stacked=show_stacked, unique=unique, save_img=save_img, l2=l2)

# @cli.command()
# def amplicon_coverage_command(samples_path: str):
#     """Get amplicon coverage"""
#     amplicon_coverage(samples_path)

# @cli.command()
# def gc_depth_command(samples_path: str):
#     """Get GC depth"""
#     gc_depth(samples_path)

if __name__ == "__main__":
    cli()
