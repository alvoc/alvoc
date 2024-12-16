import typer
from alvoc.cli.common import virus, outdir
from alvoc.core.utils import create_dir, precompute, aa, nt

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
    out = create_dir(outdir=outdir)
    genes, seq = precompute(virus, out)
    print(f"{mut} causes {aa(mut, genes, seq)[0]}")


@convert_cli.command("nt")
def convert_nt(
    virus=virus,
    mut: str = typer.Argument(
        ..., help="Nucleotide mutation in the format 'BASENPOSBASE'"
    ),
    outdir=outdir,
):
    """
    Convert nucleotide mutation to amino acid mutation for a given virus.
    """
    out = create_dir(outdir=outdir)
    genes, seq = precompute(virus, out)
    print(f"{mut} causes {nt(mut, genes, seq)}")
