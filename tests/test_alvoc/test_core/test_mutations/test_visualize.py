from pathlib import Path
from unittest.mock import patch

import pytest
from matplotlib import pyplot as plt

from alvoc.core.mutations.visualize import plot_mutations


@pytest.fixture
def mutation_data():
    sample_results = [
        {"mutation1": (10, 90), "mutation2": (20, 80)},
        {"mutation1": (15, 85), "mutation2": (25, 75)},
    ]
    sample_names = ["Sample1", "Sample2"]
    name = "test"
    min_depth = 20
    outdir = Path(".")
    return sample_results, sample_names, min_depth, name, outdir


@pytest.fixture()
def set_noninteractive_backend():
    # Set to a non-interactive backend for testing
    plt.ion()
    yield
    # Restore the original backend and interactive mode
    plt.ion()


def test_plot_generates_correctly(mutation_data, set_noninteractive_backend):
    """
    Test that the plot generates without error and handles the file path correctly.
    Ensure that no files are actually written during the test.
    """
    sample_results, sample_names, min_depth, name, outdir = mutation_data

    # Use `matplotlib.pyplot.savefig` as the target to mock
    with patch("matplotlib.pyplot.savefig") as mock_savefig:
        try:
            plot_mutations(sample_results, sample_names, min_depth, name, outdir)

            mock_savefig.assert_called_once()
        except Exception as e:
            pytest.fail(f"Plot generation failed with an error: {e}")


def test_fraction_calculation(mutation_data):
    """
    Test that mutation fractions are calculated correctly based on the input data.
    """
    sample_results, sample_names, min_depth, name, outdir = mutation_data
    expected_fractions = [
        [0.1, 0.15],  # Fractions for mutation1 across both samples
        [0.2, 0.25],  # Fractions for mutation2 across both samples
    ]
    calculated_fractions = plot_mutations(
        sample_results, sample_names, min_depth, name, outdir, return_fractions=True
    )
    assert (
        calculated_fractions == expected_fractions
    ), "Mutation fractions are not calculated correctly"
