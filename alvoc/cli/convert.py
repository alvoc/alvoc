import typer
from alvoc.cli.common import virus, outdir
from alvoc.core.utils.precompute import precompute
from alvoc.core.utils.convert import aa, nt

convert_cli = typer.Typer(no_args_is_help=True, help="Tools to convert mutations")


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
