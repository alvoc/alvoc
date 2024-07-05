from pathlib import Path
import pytest
from unittest.mock import patch, MagicMock, mock_open
from alvoc.core.mutations.prepare import prepare, process_reference, download_virus_data


def test_prepare_no_input(caplog):
    with pytest.raises(ValueError):
        prepare()
    # Check that the error log is as expected
    assert "Either 'tax_id' or 'genbank_file' must be provided." in caplog.text
    assert any("Either 'tax_id' or 'genbank_file' must be provided." in record.message for record in caplog.records)
    assert any(record.levelname == 'ERROR' for record in caplog.records)

# def test_prepare_with_genbank_file():
#     with patch('pathlib.Path.mkdir'), \
#          patch('pathlib.Path.is_dir', return_value=True):
#         assert prepare(genbank_file="/fakepath/gene_data.gb", outdir="/fakepath") is True

# def test_prepare_with_tax_id(caplog):
#     # Mock Path and file operations
#     with patch('pathlib.Path.mkdir'), \
#          patch('pathlib.Path.is_dir', return_value=True), \
#          patch('Bio.Entrez.esearch', return_value=MagicMock(read=lambda: {'IdList': ['12345']})), \
#          patch('Bio.Entrez.read', return_value={'IdList': ['67890']}), \
#          patch('Bio.Entrez.efetch', return_value=MagicMock(read=lambda: 'fake_genbank_data')), \
#          patch('Bio.Entrez.email', new_callable=lambda: 'test@example.com'), \
#          patch('builtins.open', new_callable=mock_open, create=True) as mock_file:
#         assert prepare(tax_id="2697049", outdir="/fakepath") is True
#         mock_file().write.assert_called_with('fake_genbank_data')