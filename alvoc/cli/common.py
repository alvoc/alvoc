from pathlib import Path
import typer

# commonOptions
virus: str = typer.Argument(
    ..., help="Taxonomic ID of the virus, or path to the GenBank file"
)
outdir: Path = typer.Option(
    ".", help="Output directory for results and intermediate data"
)


