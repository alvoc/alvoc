import re

CODONS = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

NTS = 'ACGT'

def aa(mut: str, genes : dict, seq : str):
    """
    Convert amino acid mutation to nucleotide mutations.

    Args:
        mut (str): The amino acid mutation in the format 'GENE:aaPOSITIONaaNEW'.
        genes (dict): Dictionary of genes with start and end positions.
        seq (str): The nucleotide sequence.

    Returns:
        list: A list of nucleotide mutations.
    """
    # Extract the gene names
    gene = mut.split(':')[0]
    
    # Handle deletion mutations
    if gene == 'DEL':
        nt_idx, length = map(int, re.findall(r'\d+', mut))
        # Convert 1-based to 0-based index
        nt_idx -= 1
        return [f'{seq[nt_idx + i]}{nt_idx + i + 1}-' for i in range(length)]
    
    # Extra the aa position
    aa_idx = int(re.findall(r'\d+', mut)[-1])
    nt_idx = genes[gene][0] + (aa_idx - 1) * 3
    codon = seq[nt_idx:nt_idx + 3]

    # Handle deletions in AA mutations
    if mut.split(':')[1].startswith('DEL') or mut.endswith('-'):
        return [f'{seq[nt_idx + i]}{nt_idx + i + 1}-' for i in range(3)]

    # Generate nucleotide mutations
    new_acid = mut[-1]
    nt_muts = []

    for i in range(4):
        if CODONS[NTS[i] + codon[1] + codon[2]] == new_acid:
            nt_muts.append(f'{codon[0]}{nt_idx + 1}{NTS[i]}')
        if CODONS[codon[0] + NTS[i] + codon[2]] == new_acid:
            nt_muts.append(f'{codon[1]}{nt_idx + 2}{NTS[i]}')
        if CODONS[codon[0] + codon[1] + NTS[i]] == new_acid:
            nt_muts.append(f'{codon[2]}{nt_idx + 3}{NTS[i]}')

    return nt_muts

def nt(mut: str, genes : dict, seq : str):
    """
    Convert nucleotide mutation to amino acid mutation.

    Args:
        mut (str): The nucleotide mutation in the format 'BASENPOSBASE'.
        genes (dict): Dictionary of genes with start and end positions.
        seq (str): The nucleotide sequence.

    Returns:
        str: The amino acid mutation.
    """             
    # Extract the mutation details
    _, new_base = mut[0], mut[-1]
    nt_idx = int(re.findall(r'\d+', mut)[0]) - 1  # Convert 1-based to 0-based index

    # Find the corresponding gene to the mutation
    for gene, (start, end) in genes.items():
        if start <= nt_idx < end:
            nt_offset = (nt_idx - start) % 3
            aa_idx = (nt_idx - start) // 3 + 1

            # Handle deletion mutations
            if new_base == '-':
                return f'{gene}:DEL{aa_idx}'

            # Generate new codon 
            codon_start = start + (aa_idx - 1) * 3
            codon = list(seq[codon_start:codon_start + 3])
            acid = CODONS[''.join(codon)]
            codon[nt_offset] = new_base
            new_acid = CODONS[''.join(codon)]
            
            return f'{gene}:{acid}{aa_idx}{new_acid}'