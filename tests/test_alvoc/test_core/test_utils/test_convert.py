import pytest
from alvoc.core.utils.convert import aa, nt

@pytest.fixture
def sample_data():
    genes = {
        'GENE1': (0, 9),    # GENE1 spans from nucleotide 0 to 8 (3 codons)
        'GENE2': (9, 18)    # GENE2 spans from nucleotide 9 to 17 (3 codons)
    }
    seq = 'ATGACCCTGAAGGCTAAGT'
    return genes, seq

def test_aa_mutation_conversion(sample_data):
    genes, seq = sample_data
    
    # Test amino acid mutation to nucleotide mutation
    mut = 'GENE1:M2T'
    expected_nt_muts = ['A4A', 'C5C', 'C6A', 'C6C', 'C6G', 'C6T']
    result = aa(mut, genes, seq)
    assert sorted(result) == sorted(expected_nt_muts), f"Unexpected result: {result}"

    # Test deletion mutation
    mut = 'DEL:4-3'
    expected_nt_muts = ['A4-', 'C5-', 'C6-']
    result = aa(mut, genes, seq)
    assert sorted(result) == sorted(expected_nt_muts), f"Unexpected result: {result}"

def test_nt_mutation_conversion(sample_data):
    genes, seq = sample_data
    
    # Test nucleotide mutation to amino acid mutation
    mut = 'A4G'
    # Position 4 in the sequence corresponds to the second codon, changing A to G
    # ACC -> GCC results in T -> A
    expected_aa_mut = 'GENE1:T2A'
    result = nt(mut, genes, seq)
    assert result == expected_aa_mut, f"Unexpected result: {result}"

    # Test deletion mutation
    mut = 'A4-'
    expected_aa_mut = 'GENE1:DEL2'
    result = nt(mut, genes, seq)
    assert result == expected_aa_mut, f"Unexpected result: {result}"
