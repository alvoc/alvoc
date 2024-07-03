import typer

from alvoc.core.amplicons import amplicon_coverage, gc_depth
from alvoc.core.mutations import aa, nt, find_mutants
from alvoc.core.lineages import find_lineages
from alvoc.core.precompute import precompute

cli = typer.Typer(help="Identify frequencies of concerning mutations from aligned reads")

@cli.callback()
def callback():
    # This function will run before any command.
    # Run precompute functions here (specific to virus of interest)
    pass 

@cli.command()
def aa_command(mut: str):
    """Process amino acid mutations"""
    aa(mut)

@cli.command()
def nt_command(mut: str):
    """Process nucleotide mutations"""
    nt(mut)

@cli.command()
def find_mutants_command(
    samples_path: str,
    mutations_path: str = typer.Option(None, "--mutations-path", "-m", help='Path to mutations'),
    min_depth: int = typer.Option(40, "--min-depth", "-d", help='Minimum depth'),
    save_img: bool = typer.Option(False, "--save-img", "-s", help='Save image'),
    csv: bool = typer.Option(False, "--csv", "-c", help='Save as CSV')
):
    """Find mutants in samples"""
    find_mutants(samples_path, mutations_path, min_depth, save_img, csv)

@cli.command()
def find_lineages_command(
    samples_path: str,
    lineages_path: str = typer.Option(None, "--lineages-path", "-l", help='Path to lineages'),
    ts: bool = typer.Option(False, "--ts", "-t", help='Time series'),
    csv: bool = typer.Option(False, "--csv", "-c", help='Save as CSV'),
    min_depth: int = typer.Option(40, "--min-depth", "-d", help='Minimum depth'),
    show_stacked: bool = typer.Option(False, "--show-stacked", "-s", help='Show stacked'),
    unique: bool = typer.Option(False, "--unique", "-u", help='Unique'),
    save_img: bool = typer.Option(False, "--save-img", "-s", help='Save image'),
    l2: bool = typer.Option(False, "--l2", "-l2", help='Level 2')
):
    """Find lineages in samples"""
    find_lineages(samples_path, lineages_path=lineages_path, ts=ts, csv=csv, min_depth=min_depth, show_stacked=show_stacked, unique=unique, save_img=save_img, l2=l2)

@cli.command()
def amplicon_coverage_command(samples_path: str):
    """Get amplicon coverage"""
    amplicon_coverage(samples_path)

@cli.command()
def gc_depth_command(samples_path: str):
    """Get GC depth"""
    gc_depth(samples_path)

if __name__ == "__main__":
    cli()
