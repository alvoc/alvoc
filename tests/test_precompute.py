import pytest
from unittest.mock import patch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from alvoc.core.utils.precompute import precompute, extract_gene_info


@pytest.fixture
def mock_seq_record():
    """Fixture to create a mock SeqRecord."""
    features = [
        SeqFeature(
            FeatureLocation(start=0, end=10),
            type="gene",
            qualifiers={"gene": ["gene1"]},
        ),
        SeqFeature(
            FeatureLocation(start=15, end=25),
            type="gene",
            qualifiers={"gene": ["gene2"]},
        ),
    ]
    return SeqRecord(Seq("ATGCATGCATGCATGCATGC"), id="test_id", features=features)


def test_extract_gene_info(mock_seq_record):
    """Test that gene information is extracted correctly."""
    expected = {
        "gene1": [0, 10],
        "gene2": [15, 25],
    }
    result = extract_gene_info(mock_seq_record)
    assert result == expected


@patch("alvoc.core.utils.precompute.process_reference")
@patch("alvoc.core.utils.precompute.download_virus_data")
def test_precompute(mock_download_virus_data, mock_process_reference, tmp_path):
    """Test precompute function end-to-end."""
    # Mock the download_virus_data and process_reference outputs
    mock_download_virus_data.return_value = tmp_path / "mock_file.gb"
    mock_process_reference.return_value = ({"gene1": (0, 10)}, "ATGC")

    virus = "Test Virus"
    outdir = tmp_path
    email = "test@example.com"

    result = precompute(virus, outdir, email)

    mock_download_virus_data.assert_called_once_with(virus, outdir, email)
    mock_process_reference.assert_called_once()
    assert result == ({"gene1": (0, 10)}, "ATGC")
