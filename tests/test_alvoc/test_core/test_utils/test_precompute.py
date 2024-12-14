from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from alvoc.core.utils.precompute import (
    download_virus_data,
    extract_gene_info,
    precompute,
    process_reference,
)

def test_precompute_with_genbank_file(tmp_path):
    virus = tmp_path / "dummy.gb"
    virus.write_text(">Dummy data")
    outdir = tmp_path / "output"

    with patch("alvoc.core.precompute.precompute.process_reference") as mocked_process:
        mocked_process.return_value = ({"gene": (1, 100)}, "ACGT", outdir)
        result = precompute(virus=str(virus), outdir=Path(outdir))
        assert result[0] == {"gene": (1, 100)}
        assert result[1] == "ACGT"
        assert result[2] == outdir


def test_precompute_with_tax_id(tmp_path, caplog):
    outdir = tmp_path / "output"
    email = "test@example.com"

    with patch(
        "alvoc.core.precompute.precompute.download_virus_data"
    ) as mocked_download:
        with patch(
            "alvoc.core.precompute.precompute.process_reference"
        ) as mocked_process:
            mocked_download.return_value = tmp_path / "downloaded.gb"
            mocked_process.return_value = ({"gene": (1, 100)}, "ACGT", outdir)
            result = precompute(virus="12345", outdir=Path(outdir), email=email)
            mocked_download.assert_called_once_with("12345", outdir, email)
            assert "No file could be processed." not in caplog.text
            assert result[0] == {"gene": (1, 100)}


def test_process_reference(tmp_path, caplog):
    reference_file = tmp_path / "dummy.gb"
    outdir = tmp_path / "output"
    outdir.mkdir()
    reference_file.write_text(">Dummy data")

    organism = SeqRecord(Seq("ACGT"), id="12345", features=[])
    with patch(
        "alvoc.core.precompute.precompute.SeqIO.parse", return_value=iter([organism])
    ):
        with patch(
            "alvoc.core.precompute.precompute.extract_gene_info",
            return_value={"gene": (1, 100)},
        ):
            result = process_reference(reference_file, outdir)
            assert result[0] == {"gene": (1, 100)}
            assert result[1] == "ACGT"
            assert result[2] == outdir
            # assert "Reference processing complete and data saved" in caplog.text


@pytest.mark.parametrize(
    "features, expected",
    [
        ([], {}),
        (
            [
                MagicMock(
                    type="gene",
                    qualifiers={"gene": ["gene1"]},
                    location=MagicMock(start=1, end=100),
                )
            ],
            {"gene1": [1, 100]},
        ),
    ],
)
def test_extract_gene_info(features, expected):
    organism = MagicMock(features=features)
    result = extract_gene_info(organism)
    assert result == expected


def test_download_virus_data_errors(caplog):
    with patch("alvoc.core.precompute.precompute.Entrez.efetch") as mocked_efetch:
        mocked_efetch.side_effect = Exception("Network error")
        with pytest.raises(Exception):
            download_virus_data("12345", Path("/fakepath"), "test@example.com")
        assert "Network error" in caplog.text
