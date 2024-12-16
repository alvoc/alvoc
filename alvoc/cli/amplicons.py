from pathlib import Path
import typer

from alvoc.cli.common import virus, outdir
from alvoc.core.utils import create_dir
from alvoc.core.utils.precompute import precompute
from alvoc.core.amplicons import Amplicons

amplicons_cli = typer.Typer(no_args_is_help=True, help="Tools to analyze amplicons")


@amplicons_cli.command()
def coverage(
    virus=virus,
    samples: Path = typer.Argument(
        ..., help="Path to a BAM file or CSV file listing samples."
    ),
    inserts_path: Path = typer.Argument(
        ..., help="Path to the CSV file detailing the regions (amplicons) to evaluate."
    ),
    outdir=outdir,
):
    """Get amplicon coverage"""
    out = create_dir(outdir=outdir)
    _, seq= precompute(virus, out)
    Amplicons(
        file_path=samples,
        inserts=inserts_path,
        sequence=seq,
        outdir=out,
    ).amplicon_coverage()


@amplicons_cli.command()
def gc_depth(
    virus=virus,
    samples: Path = typer.Argument(
        ..., help="Path to a BAM file or CSV file listing samples."
    ),
    inserts_path: Path = typer.Argument(
        ..., help="Path to the CSV file detailing the regions (amplicons) to evaluate."
    ),
    outdir=outdir,
):
    """Get GC depth"""
    out = create_dir(outdir=outdir)
    _, seq= precompute(virus, out)
    Amplicons(
        file_path=samples,
        inserts=inserts_path,
        sequence=seq,
        outdir=out,
    ).gc_depth()
