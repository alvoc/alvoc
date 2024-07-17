from pathlib import Path
import typer

from alvoc.core.logging import init_logger
from alvoc.core.mutations.analyze import find_mutants as find_mutants_internal
from alvoc.core.precompute import precompute
from alvoc.core.utils.convert import aa, nt
from alvoc.core.lineages import find_lineages
from alvoc.core.amplicons import Amplicons

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
    samples_path: str = typer.Argument(
        ..., help="Path to the file listing samples or a single BAM file."
    ),
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


@cli.command()
def find_lineages_command(
    samples_path: str = typer.Argument(
        ..., help="Path to the file listing samples or a single BAM file."
    ),
    tax_id: str = typer.Option(
        None, help="Taxonomic ID of the virus, required if no GenBank file is provided"
    ),
    genbank_file: str = typer.Option(
        None, help="Path to the GenBank file, required if no tax ID is provided"
    ),
    lineages_path: str = typer.Option(
        None, "--lineages-path", "-l", help="Path to lineages"
    ),
    black_list: list[str] = typer.Option(
        None, "--black-list", "-b", help="List of lineages to black list"
    ),
    ts: bool = typer.Option(False, "--ts", "-t", help="Time series"),
    min_depth: int = typer.Option(40, "--min-depth", "-d", help="Minimum depth"),
    show_stacked: bool = typer.Option(
        False, "--show-stacked", "-s", help="Show stacked"
    ),
    unique: bool = typer.Option(False, "--unique", "-u", help="Unique"),
    l2: bool = typer.Option(False, "--l2", "-l2", help="Level 2"),
    outdir: str = typer.Option(
        ".", help="Output directory for results and intermediate data"
    ),
):
    """Find lineages in samples"""
    genes, seq, out = precompute(tax_id, genbank_file, outdir)
    mut_lins = {}
    file_path = Path(samples_path)
    lineages_file = Path(lineages_path)
    find_lineages(
        file_path=file_path,
        mut_lins=mut_lins,
        genes=genes,
        seq=seq,
        outdir=out,
        black_list=black_list,
        lineages_path=lineages_file,
        ts=ts,
        min_depth=min_depth,
        show_stacked=show_stacked,
        unique=unique,
        l2=l2,
    )


@cli.command()
def amplicon_coverage(
    samples_path: str = typer.Argument(
        ..., help="Path to the file listing samples or a single BAM file."
    ),
    inserts_path: str = typer.Argument(
        ..., help="Path to the CSV file detailing the regions (amplicons) to evaluate."
    ),
    tax_id: str = typer.Option(
        None, help="Taxonomic ID of the virus, required if no GenBank file is provided"
    ),
    genbank_file: str = typer.Option(
        None, help="Path to the GenBank file, required if no tax ID is provided"
    ),
    outdir: str = typer.Option(
        ".", help="Output directory for results and intermediate data"
    ),
):
    """Get amplicon coverage"""
    _, seq, out = precompute(tax_id, genbank_file, outdir)
    Amplicons(
        file_path=Path(samples_path),
        inserts=Path(inserts_path),
        sequence=seq,
        outdir=out,
    ).amplicon_coverage()


@cli.command()
def gc_depth(
    samples_path: str = typer.Argument(
        ..., help="Path to the file listing samples or a single BAM file."
    ),
    inserts_path: str = typer.Argument(
        ..., help="Path to the CSV file detailing the regions (amplicons) to evaluate."
    ),
    tax_id: str = typer.Option(
        None, help="Taxonomic ID of the virus, required if no GenBank file is provided"
    ),
    genbank_file: str = typer.Option(
        None, help="Path to the GenBank file, required if no tax ID is provided"
    ),
    outdir: str = typer.Option(
        ".", help="Output directory for results and intermediate data"
    ),
):
    """Get GC depth"""
    _, seq, out = precompute(tax_id, genbank_file, outdir)
    Amplicons(
        file_path=Path(samples_path),
        inserts=Path(inserts_path),
        sequence=seq,
        outdir=out,
    ).gc_depth()



if __name__ == "__main__":
    cli()
