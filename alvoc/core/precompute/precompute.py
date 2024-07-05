import json
from Bio import SeqIO, Entrez
from pathlib import Path
from alvoc.core import logging

logger = logging.get_logger()


def precompute(
    tax_id: str | None = None,
    genbank_file: str | None = None,
    outdir: str = "",
    email: str = "example@example.com",
) -> tuple[dict[str, tuple[int, int]], str, Path]:
    """
    Processes a GenBank file to extract gene information and sequence, or alternatively pass in a virus taxonomic ID
    to automatically generate the necessary reference data. At least one of 'tax_id' or 'genbank_file' must be provided.

    Args:
        tax_id : Taxonomic ID of the virus. Required if 'genbank_file' is not provided.
        genbank_file : Path to the GenBank file. Required if 'tax_id' is not provided.
        email : Email for accessing Entrez api.
        outdir : Output directory for results and intermediate data. Defaults to the current directory.

    Returns:
        A tuple with gene coordinates (dictionary) and genome sequence (string).

    Raises:
        ValueError: If neither 'tax_id' nor 'genbank_file' is provided.
    """
    if not (tax_id or genbank_file):
        raise ValueError("Either 'tax_id' or 'genbank_file' must be provided.")

    outdir_path = Path(outdir)
    if not outdir_path.is_dir():
        logger.info(f"Creating directory at {outdir_path}")
        outdir_path.mkdir(parents=True, exist_ok=True)

    if genbank_file:
        reference_file = Path(genbank_file)
        return process_reference(reference_file, outdir_path)
    elif tax_id:
        file_path = download_virus_data(tax_id, outdir_path, email)
        if file_path:
            reference_file = Path(file_path)
            return process_reference(reference_file, outdir_path)
        else:
            logger.error("No file could be processed.")
            raise ValueError("No file could be processed.")
    else:
        logger.error("Either 'tax_id' or 'genbank_file' must be provided.")
        raise ValueError("Either 'tax_id' or 'genbank_file' must be provided.")


def process_reference(
    reference_file: Path, outdir_path
) -> tuple[dict[str, tuple[int, int]], str, Path]:
    """
    Processes a GenBank file to extract gene information and sequence.

    Args:
        reference_file: Path object with Genbank file.
        outdir_path : Path object with outdir directory.

    Returns:
        A tuple with gene coordinates (dictionary) and genome sequence (string).
    """
    logger.info("Processing reference")

    try:
        organism = next(SeqIO.parse(reference_file.as_posix(), "genbank"))
        gene_coordinates = extract_gene_info(organism)
        genome_sequence = str(organism.seq)

        file_path = outdir_path / "gene_data.json"

        with open(file_path, "w") as f:
            json.dump(
                {
                    "gene_coordinates": gene_coordinates,
                    "genome_sequence": genome_sequence,
                },
                f,
                indent=4,
            )
        logger.info("Reference processing complete and data saved")

        return gene_coordinates, genome_sequence, outdir_path

    except Exception as e:
        logger.error(f"An error occurred: {e}")
        return {}, "", Path(".")


def extract_gene_info(organism):
    gene_coordinates = {}
    for feature in organism.features:
        if feature.type == "gene":
            try:
                gene = feature.qualifiers["gene"][0]
                start = int(feature.location.start)
                end = int(feature.location.end)
                gene_coordinates[gene] = [start, end]
            except KeyError:
                logger.info(f"Skipping feature with no 'gene' qualifier: {feature}")
    return gene_coordinates


def download_virus_data(tax_id: str, outdir: Path, email: str):
    """
    Downloads virus gene data from GenBank based on a given taxonomic ID.

    Args:
        tax_id : Taxonomic ID of the virus.
        outdir : Path to the output directory.
        email : Email for Entrez API.
    """
    try:
        Entrez.email = email
        search_handle = Entrez.esearch(
            db="nucleotide", term=f"txid{tax_id}[Organism:exp]", retmax=1
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if (
            not isinstance(search_results, dict)
            or "IdList" not in search_results
            or not search_results["IdList"]
        ):
            raise Exception(f"No results found for tax ID: {tax_id}")

        virus_id = search_results["IdList"][0]
        fetch_handle = Entrez.efetch(
            db="nucleotide", id=virus_id, rettype="gb", retmode="text"
        )

        file_path = outdir / "gene_data.gb"
        with open(file_path, "w") as f:
            f.write(fetch_handle.read())
        fetch_handle.close()
        logger.info(f"Downloaded data for tax ID {tax_id} and saved to {file_path}")
        return file_path.as_posix()
    except Exception as e:
        logger.error(f"An error occurred while downloading data: {e}")
        raise