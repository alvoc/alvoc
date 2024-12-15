import pytest
from pathlib import Path
from unittest.mock import patch
from alvoc.core.amplicons.analyze import Amplicons
from pandas import DataFrame


def test_initialization_correct_columns():
    """
    Test initialization of Amplicons class with correct CSV formats.
    """
    with patch(
        "alvoc.core.amplicons.analyze.pd.read_csv",
        return_value=DataFrame(columns=["chrom", "chromStart", "chromEnd"]),
    ):
        # Initialization should succeed without exceptions
        analyzer = Amplicons(
            Path("path/to/file.bam"),
            Path("path/to/inserts.csv"),
            "AGCT" * 1000,
            Path("/output"),
        )
        assert analyzer is not None


def test_initialization_incorrect_columns():
    """
    Test initialization of Amplicons class with incorrect CSV formats.
    """
    with patch(
        "alvoc.core.amplicons.analyze.pd.read_csv",
        return_value=DataFrame(columns=["wrong1", "wrong2"]),
    ):
        # Expect ValueError due to missing required columns
        with pytest.raises(
            ValueError,
            match="CSV file must contain 'chrom', 'chromStart', 'chromEnd' columns.",
        ):
            Amplicons(
                Path("path/to/file.bam"),
                Path("path/to/inserts.csv"),
                "AGCT" * 1000,
                Path("/output"),
            )


def test_calculate_gc_content():
    """
    Test the static method calculate_gc_content.
    """
    gc_content = Amplicons.calculate_gc_content("GCGC")
    assert gc_content == 1.0  # 100% GC content

    gc_content = Amplicons.calculate_gc_content("ATAT")
    assert gc_content == 0.0  # 0% GC content
