import typer
from pathlib import Path

from alvoc.cli.common import outdir
from alvoc.core.constellations import (
    make_constellations,
    NextstrainSource,
    MSADataSource,
)
from alvoc.core.utils import create_dir

constellations_cli = typer.Typer(
    no_args_is_help=True, help="Tools to make constellations"
)


@constellations_cli.command("nextstrain")
def nextstrain(
    tree_url: str = typer.Argument(..., help="Nextstrain phylogeny tree dataset JSON URL"),
    proportion_threshold: int = typer.Option(
        0.9,
        "--proportion_threshold",
        "-pt",
        help="Minimum proportion of nodes in a clade required to include a mutation",
    ),
    outdir=outdir,
):
    """
    Generates constellations using the provided nextstrain phylogeny dataset.
    """
    out = create_dir(outdir=outdir)
    src = NextstrainSource()
    make_constellations(
        source=src, source_path=tree_url, outdir=out, threshold=proportion_threshold
    )


@constellations_cli.command("msa")
def msa(
    fasta: Path = typer.Argument(..., exists=True, help="MSA FASTA file"),
    proportion_threshold: int = typer.Option(
        0.9,
        "--proportion_threshold",
        "-pt",
        help="Minimum proportion of nodes in a clade required to include a mutation",
    ),
    outdir=outdir,
):
    """
    Generates constellations using a custom MSA FASTA file.
    """
    out = create_dir(outdir=outdir)
    src = MSADataSource()
    make_constellations(
        source=src, source_path=str(fasta), outdir=out, threshold=proportion_threshold
    )


if __name__ == "__main__":
    constellations_cli()
