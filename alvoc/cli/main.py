import typer
import json
from pathlib import Path

from alvoc.core.utils.precompute import precompute
from alvoc.cli.common import virus, outdir, with_spinner
from alvoc.core.mutations.analyze import find_mutants as find_mutants_internal

from alvoc.core.lineages import find_lineages as fl
from alvoc.core.utils.logging import init_logger

from alvoc.cli.amplicons import amplicons_cli
from alvoc.cli.convert import convert_cli
from importlib.metadata import version as get_version

cli = typer.Typer(
    no_args_is_help=True,
    help="Identify frequencies of concerning mutations from aligned reads",
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

cli.add_typer(amplicons_cli, name="amplicons")
cli.add_typer(convert_cli, name="convert")


def version_callback(value: bool):
    if value:
        typer.echo(f"Alvoc: {get_version("alvoc")}")
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
        ..., help="Path to a BAM file or TXT file listing samples."
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


if __name__ == "__main__":
    cli()
