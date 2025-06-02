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
    proportion_threshold: float = typer.Option(
        0.9,
        "--proportion_threshold",
        "-pt",
        help="Minimum proportion of nodes in a clade required to include a mutation",
    ),
    clade_delim: str = typer.Option(
        None,
        "--clade-delim",
        "-cd",
        help="Delimiter character/string to split each FASTA header on",
    ),
    clade_field: int = typer.Option(
        None,
        "--clade-field",
        "-cf",
        help="Zero-based field index, after splitting on --clade-delim, that contains the clade name",
    ),
    outdir: Path = typer.Option(
        Path("."), "--outdir", "-o", help="Output directory for constellations"
    ),
):
    """
    Generates constellations using a custom MSA FASTA file.
    """
    # 1) Create output directory (like you already do)
    out = outdir
    out.mkdir(parents=True, exist_ok=True)

    # 2) Pass the user’s delimiter + field into MSADataSource
    src = MSADataSource(clade_delim=clade_delim, clade_field=clade_field)

    # 3) Call into the generic pipeline
    make_constellations(
        source=src,
        source_path=str(fasta),
        outdir=out,
        threshold=proportion_threshold,
    )
 