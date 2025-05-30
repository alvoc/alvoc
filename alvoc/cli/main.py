import typer
from pathlib import Path

from alvoc.cli.common import virus, outdir, with_spinner
from alvoc.core.amplicons.main import calculate_amplicon_metrics
from alvoc.core.utils import create_dir
from alvoc.core.utils.precompute import precompute
from alvoc.core.variants.mutations import find_mutants as fm

from alvoc.core.variants.lineages import find_lineages as fl
from alvoc.core.utils.logging import init_logger

from alvoc.cli.convert import convert_cli
from alvoc.cli.constellations import constellations_cli
from importlib.metadata import version as get_version
from typer.main import get_command

cli = typer.Typer(
    no_args_is_help=True,
    help="Abundance learning for variants of concern",
)

# Inject spinner into all commands
original_command = cli.command


def command_with_spinner(*args, **kwargs):
    def decorator(func):
        # Apply the original command decorator first
        decorated_command = original_command(*args, **kwargs)
        # Then apply the spinner decorator
        return decorated_command(with_spinner(func))

    return decorator


cli.command = command_with_spinner

cli.add_typer(convert_cli, name="convert")
cli.add_typer(constellations_cli, name="constellations")


def version_callback(value: bool):
    if value:
        typer.echo(f"Alvoc: {get_version('alvoc')}")
        raise typer.Exit()


@cli.callback()
def callback(
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        help="Current version of Alvoc",
        is_eager=True,
    ),
):
    # Set default logging to INFO
    init_logger(20)
    pass


@cli.command()
def find_lineages(
    virus=virus,
    samples: Path = typer.Argument(
        ..., help="Path to a BAM file or CSV file listing samples."
    ),
    constellations: Path = typer.Argument(
        ..., help="Path to a JSON file containing mutation lineage constellations."
    ),
    white_list: Path = typer.Option(
        None,
        "--white-list",
        "-wl",
        help="Path to a TXT file containing lineages to inclue.",
    ),
    black_list: Path = typer.Option(
        None,
        "--black-list",
        "-bl",
        help="Path to a TXT file containing lineages to exclude.",
    ),
    min_depth: int = typer.Option(
        40, "--min-depth", "-d", help="Minimum depth for a mutation to be considered."
    ),
    unique: bool = typer.Option(
        False, "--unique", "-u", help="Whether to consider unique mutations only."
    ),
    l2: bool = typer.Option(
        False,
        "--l2",
        "-l2",
        help="Whether to use a secondary method for regression analysis.",
    ),
    outdir=outdir,
):
    """Find lineages in samples"""
    fl(
        virus=virus,
        samples=samples,
        constellations=constellations,
        outdir=outdir,
        white_list=white_list,
        black_list=black_list,
        min_depth=min_depth,
        unique=unique,
        l2=l2,
    )


@cli.command()
def find_mutants(
    virus=virus,
    samples: Path = typer.Argument(
        ..., help="Path to a BAM file or CSV file listing samples."
    ),
    constellations: Path = typer.Argument(
        ..., help="Path to a JSON file containing mutation lineage constellations."
    ),
    mutations: Path = typer.Argument(..., help="Path to mutations"),
    min_depth: int = typer.Option(10, "--min-depth", "-d", help="Minimum depth"),
    outdir=outdir,
):
    """
    Find mutations in sequencing data, either from BAM files or a sample list.
    """
    fm(
        virus=virus,
        samples=samples,
        constellations=constellations,
        mutations_path=mutations,
        min_depth=min_depth,
        outdir=outdir,
    )


@cli.command()
def extract_gene_data(virus=virus, outdir=outdir):
    """
    Extracts gene coordinates and the genome sequence from a GenBank file or generates them using an Entrez API search term.
    """
    out = create_dir(outdir=outdir)
    precompute(virus=virus, outdir=out)


@cli.command()
def amplicons(
    virus=virus,
    samples: Path = typer.Argument(
        ..., help="Path to a BAM file or CSV file listing samples."
    ),
    inserts_path: Path = typer.Argument(
        ..., help="Path to the BED file detailing the regions (amplicons) to evaluate."
    ),
    outdir=outdir,
):
    """Get amplicon metrics such as coverage, gc_content and visualizations"""
    calculate_amplicon_metrics(virus, samples, inserts_path, outdir)


click_cli = get_command(cli)


if __name__ == "__main__":
    cli()
