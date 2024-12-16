import json
from pathlib import Path

from Bio import Entrez, SeqIO

from alvoc.core.utils.logging import get_logger

logger = get_logger()


def precompute(
    virus: str,
    outdir: Path,
    email: str = "example@example.com",
) -> tuple[dict[str, tuple[int, int]], str]:
    """
    Processes a GenBank file to extract gene information and sequence, or alternatively pass in an NCBI query to automatically generate the necessary reference data.

    Args:
        virus : NCBI query for virus or Path to the GenBank file
        email : Email for accessing Entrez api.
        outdir : Output directory for results and intermediate data. Defaults to the current directory.

    Returns:
        A tuple with gene coordinates (dictionary) and genome sequence (string).

    """
    reference_file = Path(virus)
    if reference_file.exists():
        return process_reference(reference_file, outdir)
    else:
        file_path = download_virus_data(virus, outdir, email)
        if file_path:
            reference_file = Path(file_path)
            return process_reference(reference_file, outdir)
        else:
            logger.error("No file could be processed.")
            raise ValueError("No file could be processed.")


def process_reference(
    reference_file: Path, outdir_path: Path
) -> tuple[dict[str, tuple[int, int]], str]:
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

        return gene_coordinates, genome_sequence

    except Exception as e:
        logger.error(f"An error occurred: {e}")
        return {}, ""


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


def download_virus_data(query: str, outdir: Path, email: str):
    """
    Downloads virus gene data from GenBank based on a given taxonomic ID.

    Args:
        query : Search term for the virus.
        outdir : Path to the output directory.
        email : Email for Entrez API.
    """
    try:
        Entrez.email = email
        search_handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if (
            not isinstance(search_results, dict)
            or "IdList" not in search_results
            or not search_results["IdList"]
        ):
            raise Exception(f"No results found for: {query}")

        virus_id = search_results["IdList"][0]
        fetch_handle = Entrez.efetch(
            db="nucleotide", id=virus_id, rettype="gb", retmode="text"
        )

        file_path = outdir / "gene_data.gb"
        with open(file_path, "w") as f:
            f.write(fetch_handle.read())
        fetch_handle.close()
        logger.info(f"Downloaded data for query and saved to {file_path}")
        return file_path.as_posix()
    except Exception as e:
        logger.error(f"An error occurred while downloading data: {e}")
        raise
