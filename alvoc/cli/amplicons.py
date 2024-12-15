from pathlib import Path
import typer

from alvoc.cli.common import virus, outdir
from alvoc.core.utils.precompute import precompute
from alvoc.core.amplicons import Amplicons

amplicons_cli = typer.Typer(
    no_args_is_help=True, help="Tools to analyze amplicons"
)

@amplicons_cli.command()
def coverage(
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


@amplicons_cli.command()
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

