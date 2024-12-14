import json
from pathlib import Path
import typer

from alvoc.core.utils.logging import init_logger
from alvoc.core.mutations.analyze import find_mutants as find_mutants_internal
from alvoc.core.utils.precompute import precompute
from alvoc.core.utils.convert import aa, nt
from alvoc.core.lineages import find_lineages as fl
from alvoc.core.amplicons import Amplicons

cli = typer.Typer(
    no_args_is_help=True,
    help="Identify frequencies of concerning mutations from aligned reads",
)


@cli.callback()
def callback():
    init_logger()
    pass


# commonOptions
virus: str = typer.Argument(
    ..., help="Taxonomic ID of the virus, or path to the GenBank file"
)
outdir: Path = typer.Option(
    ".", help="Output directory for results and intermediate data"
)

convert_cli = typer.Typer(
    no_args_is_help=True, help="Tools to convert mutations"
)
cli.add_typer(convert_cli, name="convert")


@convert_cli.command("aa")
def convert_aa(
    virus=virus,
    mut: str = typer.Argument(
        ..., help="Amino acid mutation in the format 'GENE:aaPOSITIONaaNEW'"
    ),
    outdir=outdir,
):
    """
    Convert amino acid mutation to nucleotide mutations for a given virus.
    """
    genes, seq, _ = precompute(virus, outdir)
    return aa(mut, genes, seq)


@convert_cli.command("nt")
def convert_nt(
    virus=virus,
    mut: str = typer.Option(
        ..., help="Nucleotide mutation in the format 'BASENPOSBASE'"
    ),
    outdir=outdir,
):
    """
    Convert nucleotide mutation to amino acid mutation for a given virus.
    """
    genes, seq, _ = precompute(virus, outdir)
    return nt(mut, genes, seq)


@cli.command()
def find_mutants(
    virus=virus,
    samples_path: Path = typer.Argument(
        ..., help="Path to the file listing samples or a single BAM file."
    ),
    mut_lin_path: Path = typer.Argument(
        ..., help="Path to the json file containing mutation lineages."
    ),
    mutations_path: Path = typer.Argument(..., help="Path to mutations"),
    min_depth: int = typer.Option(40, "--min-depth", "-d", help="Minimum depth"),
    outdir=outdir,
):
    """
    Find mutations in sequencing data, either from BAM files or a sample list.
    """
    genes, seq, out = precompute(virus, outdir)
    mut_lins = {}
    with open(mut_lin_path, "r") as file:
        mut_lins = json.load(file)
    find_mutants_internal(
        samples_path, mutations_path, min_depth, mut_lins, genes, seq, out
    )


@cli.command()
def find_lineages(
    virus=virus,
    samples_path: Path = typer.Argument(
        ..., help="Path to the file listing samples or a single BAM file."
    ),
    mut_lin_path: Path = typer.Argument(
        ..., help="Path to the json file containing mutation lineages."
    ),
    lineages_path: Path = typer.Option(
        None, "--lineages-path", "-l", help="Path to lineages"
    ),
    black_list: list[str] = typer.Option(
        [], "--black-list", "-b", help="List of lineages to black list"
    ),
    ts: bool = typer.Option(False, "--ts", "-t", help="Time series"),
    min_depth: int = typer.Option(40, "--min-depth", "-d", help="Minimum depth"),
    show_stacked: bool = typer.Option(
        False, "--show-stacked", "-s", help="Show stacked"
    ),
    unique: bool = typer.Option(False, "--unique", "-u", help="Unique"),
    l2: bool = typer.Option(False, "--l2", "-l2", help="Level 2"),
    outdir=outdir,
):
    """Find lineages in samples"""
    genes, seq, out = precompute(virus, outdir)
    fl(
        file_path=samples_path,
        mut_lin_path=mut_lin_path,
        genes=genes,
        seq=seq,
        outdir=out,
        black_list=black_list,
        lineages_path=lineages_path,
        ts=ts,
        min_depth=min_depth,
        show_stacked=show_stacked,
        unique=unique,
        l2=l2,
    )


@cli.command()
def amplicon_coverage(
    virus=virus,
    samples_path: str = typer.Argument(
        ..., help="Path to the file listing samples or a single BAM file."
    ),
    inserts_path: str = typer.Argument(
        ..., help="Path to the CSV file detailing the regions (amplicons) to evaluate."
    ),
    outdir=outdir,
):
    """Get amplicon coverage"""
    _, seq, out = precompute(virus, outdir)
    Amplicons(
        file_path=Path(samples_path),
        inserts=Path(inserts_path),
        sequence=seq,
        outdir=out,
    ).amplicon_coverage()


@cli.command()
def gc_depth(
    virus=virus,
    samples_path: str = typer.Argument(
        ..., help="Path to the file listing samples or a single BAM file."
    ),
    inserts_path: str = typer.Argument(
        ..., help="Path to the CSV file detailing the regions (amplicons) to evaluate."
    ),
    outdir=outdir,
):
    """Get GC depth"""
    _, seq, out = precompute(virus, outdir)
    Amplicons(
        file_path=Path(samples_path),
        inserts=Path(inserts_path),
        sequence=seq,
        outdir=out,
    ).gc_depth()


if __name__ == "__main__":
    cli()
