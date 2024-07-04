# import pytest
# from unittest.mock import patch, mock_open
# from alvoc.core.mutations.analyze import find_mutants

# # Fixture for mutation data and sample data
# @pytest.fixture
# def mutation_data():
#     return {
#         'mutations.txt': "mut1\nmut2",
#         'sample1.bam': ['mut1', 'mut2'],
#     }

# @pytest.fixture
# def sample_file_data():
#     return "sample1.bam\tSample1\nsample2.bam\tSample2"

# # Mock for reading mutation and sample files
# @pytest.fixture
# def mock_file_open(mutation_data, sample_file_data):
#     m = mock_open()
#     m.side_effect = [mock_open(read_data=mutation_data['mutations.txt']).return_value,
#                      mock_open(read_data=sample_file_data).return_value]
#     return m

# # Test for the find_mutants function
# def test_find_mutants(mutation_data, mock_file_open):
#     with patch("builtins.open", mock_file_open), \
#          patch("pysam.AlignmentFile") as mock_samfile, \
#          patch("alvoc.core.mutations.analyze.find_mutants_in_bam") as mock_find_mutants_in_bam:

#         # Set up mocks
#         mock_find_mutants_in_bam.return_value = {'mut1': [10, 90], 'mut2': [5, 95]}

#         # Execute the function
#         mut_lins = {'lineage1': {'mut1': 1, 'mut2': 0}}
#         find_mutants('sample_list.txt', 'mutations.txt', 10, False, False, mut_lins)

#         # Assertions to check if the function behaves as expected
#         mock_find_mutants_in_bam.assert_called_with('sample1.bam', ['mut1', 'mut2'])
#         mock_samfile.assert_called_with('sample1.bam', 'rb')
#         # Add more assertions here based on your actual function requirements and expected behaviors

