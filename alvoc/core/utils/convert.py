import re

CODONS = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "_",
    "TAG": "_",
    "TGC": "C",
    "TGT": "C",
    "TGA": "_",
    "TGG": "W",
}

NTS = "ACGT"


def aa(mut: str, genes: dict, seq: str) -> list[str]:
    """
    Convert an amino acid mutation descriptor to corresponding nucleotide mutations.

    Args:
        mut: Amino acid mutation in the format 'GENE:aaPOSITIONaaNEW'.
        genes: Dictionary mapping gene names to their start and end positions in the sequence.
        seq: Complete nucleotide sequence as a string.

    Returns:
        A list of nucleotide mutations derived from the specified amino acid mutation.
    """

    # Extract the gene names
    gene = mut.split(":")[0]

    # Handle deletion mutations
    if gene == "DEL":
        nt_idx, length = map(int, re.findall(r"\d+", mut))
        # Convert 1-based to 0-based index
        nt_idx -= 1
        return [f"{seq[nt_idx + i]}{nt_idx + i + 1}-" for i in range(length)]

    # Extra the aa position
    aa_idx = int(re.findall(r"\d+", mut)[-1])
    nt_idx = genes[gene][0] + (aa_idx - 1) * 3
    codon = seq[nt_idx : nt_idx + 3]

    # Handle deletions in AA mutations
    if mut.split(":")[1].startswith("DEL") or mut.endswith("-"):
        return [f"{seq[nt_idx + i]}{nt_idx + i + 1}-" for i in range(3)]

    # Generate nucleotide mutations
    new_acid = mut[-1]
    nt_muts = []

    for i in range(4):
        if CODONS[NTS[i] + codon[1] + codon[2]] == new_acid:
            nt_muts.append(f"{codon[0]}{nt_idx + 1}{NTS[i]}")
        if CODONS[codon[0] + NTS[i] + codon[2]] == new_acid:
            nt_muts.append(f"{codon[1]}{nt_idx + 2}{NTS[i]}")
        if CODONS[codon[0] + codon[1] + NTS[i]] == new_acid:
            nt_muts.append(f"{codon[2]}{nt_idx + 3}{NTS[i]}")

    return nt_muts


def nt(mut: str, genes: dict, seq: str) -> str:
    """
    Convert nucleotide mutation to amino acid mutation.

    Args:
        mut : The nucleotide mutation in the format 'BASENPOSBASE'.
        genes : Dictionary of genes with start and end positions.
        seq : The nucleotide sequence.

    Returns:
        The amino acid mutation derived from the specified nucleotide mutation.

    Raises:
        ValueError: If the mutation index is outside of any gene regions or if input format is invalid.
    """
    if not re.match(r"^[ATCG]-?\d+[ATCG-]$", mut):
        raise ValueError("Invalid mutation format.")

    _, new_base = mut[0], mut[-1]
    nt_idx = int(re.findall(r"\d+", mut)[0]) - 1  # Convert 1-based to 0-based index

    # Find the corresponding gene to the mutation
    for gene, (start, end) in genes.items():
        if start <= nt_idx < end:
            nt_offset = (nt_idx - start) % 3
            aa_idx = (nt_idx - start) // 3 + 1

            # Handle deletion mutations
            if new_base == "-":
                return f"{gene}:DEL{aa_idx}"

            # Generate new codon
            codon_start = start + (aa_idx - 1) * 3
            codon = list(seq[codon_start : codon_start + 3])
            codon[nt_offset] = new_base
            new_acid = CODONS.get(
                "".join(codon), "?"
            )  # '?' as a placeholder for unknown codon translation

            acid = CODONS.get("".join(seq[codon_start : codon_start + 3]), "?")

            return f"{gene}:{acid}{aa_idx}{new_acid}"

    # Handle case where the mutation index does not correspond to any gene
    raise ValueError("Mutation index is outside of the gene regions.")
