import pytest
from alvoc.core.mutations.visualize import plot_mutations
import matplotlib
import os


@pytest.fixture
def mutation_data():
    sample_results = [
        {"mutation1": (10, 90), "mutation2": (20, 80)},
        {"mutation1": (15, 85), "mutation2": (25, 75)},
    ]
    sample_names = ["Sample1", "Sample2"]
    min_depth = 20
    img_path = (
        None  # Set to None to show plot or 'path/to/your/test/output.png' to save
    )
    return sample_results, sample_names, min_depth, img_path


def test_plot_generates_correctly(mutation_data):
    """
    Test that the plot generates without error and handles the file path correctly.
    """
    # Ensure plot is not interactive to prevent issues with CI
    matplotlib.use("Agg")
    sample_results, sample_names, min_depth, img_path = mutation_data
    try:
        plot_mutations(sample_results, sample_names, min_depth, img_path)
        assert True  # If no error, the plot generation is assumed successful
    except Exception as e:
        pytest.fail(f"Plot generation failed with an error: {e}")


def test_plot_saves_file(mutation_data):
    """
    Test that the plot saves to a file if img_path is provided.
    """
    sample_results, sample_names, min_depth, _ = mutation_data
    test_path = "test_output.png"
    plot_mutations(sample_results, sample_names, min_depth, test_path)
    assert os.path.isfile(test_path), "File was not created where expected"
    os.remove(test_path)  # Cleanup after test


def test_fraction_calculation(mutation_data):
    """
    Test that mutation fractions are calculated correctly based on the input data.
    """
    sample_results, sample_names, min_depth, _ = mutation_data
    expected_fractions = [
        [0.1, 0.15],  # Fractions for mutation1 across both samples
        [0.2, 0.25],  # Fractions for mutation2 across both samples
    ]
    calculated_fractions = plot_mutations(
        sample_results, sample_names, min_depth, None, return_fractions=True
    )
    assert (
        calculated_fractions == expected_fractions
    ), "Mutation fractions are not calculated correctly"
