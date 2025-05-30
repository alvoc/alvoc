# msa.py
from typing import Iterator, Tuple, List
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from alvoc.core.constellations.core import DataSource


class MSADataSource(DataSource):
    """
    Multiple‐Sequence‐Alignment data source.
    Uses the first sequence in the alignment as reference,
    and yields one record per query sequence.
    """

    def fetch(self, source: str) -> MultipleSeqAlignment:
        """
        Load a FASTA MSA from `source` (file path or URL).
        """
        # (if you need to support URLs, you could download to a temp file here)
        return AlignIO.read(source, "fasta")

    def records(
        self, raw_data: MultipleSeqAlignment
    ) -> Iterator[Tuple[str, List[str]]]:
        """
        Compare each sequence to the reference (alignment[0]),
        build mutations like "A100T", and yield (sequence_id, [mutations]).
        """
        # build a map from alignment column → ungapped reference position
        ref_seq = str(raw_data[0].seq)
        ref_map: List[int] = []
        counter = 0
        for base in ref_seq:
            if base != "-":
                counter += 1
                ref_map.append(counter)
            else:
                ref_map.append(None)

        # for each non-reference record:
        for record in raw_data[1:]:
            seq_id = record.id
            alt_seq = str(record.seq)
            muts: List[str] = []

            for idx, (r_base, a_base) in enumerate(zip(ref_seq, alt_seq)):
                # skip gaps
                if r_base == "-" or a_base == "-":
                    continue
                # record a mutation
                if r_base != a_base:
                    pos = ref_map[idx]  # guaranteed non‐None here
                    muts.append(f"{r_base}{pos}{a_base}")

            # yield one profile per sequence
            yield seq_id, muts
