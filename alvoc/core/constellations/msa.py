# msa.py
from typing import Iterator, Tuple, List
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from alvoc.core.constellations.core import DataSource

class MSADataSource(DataSource):
    """
    MSADataSource now takes:
      • clade_delim  (e.g. "|" or "," or ";")
      • clade_field  (zero-based index of which split-piece is the clade)
    If both are provided, we do header.split(clade_delim)[clade_field].
    Otherwise, we leave `clade_key = full_header`.
    """

    def __init__(self, clade_delim: str = None, clade_field: int = None):
        self.clade_delim = clade_delim
        self.clade_field = clade_field

    def fetch(self, source: str) -> MultipleSeqAlignment:
        return AlignIO.read(source, "fasta")

    def records(
        self, raw_data: MultipleSeqAlignment
    ) -> Iterator[Tuple[str, List[str]]]:
        # 1) Build a map from alignment-column → reference-position:
        ref_seq = str(raw_data[0].seq)
        ref_map: List[int] = []
        counter = 0
        for base in ref_seq:
            if base != "-":
                counter += 1
                ref_map.append(counter)
            else:
                ref_map.append(None)

        # 2) Loop over each non-reference record
        for record in raw_data[1:]:
            header = record.id
            # If both delim+field are set, split; otherwise use the full header
            if self.clade_delim is not None and self.clade_field is not None:
                parts = header.split(self.clade_delim)
                # Guard against out-of-range indices:
                if self.clade_field < len(parts):
                    clade_key = parts[self.clade_field]
                else:
                    # fallback if the user’s choices are invalid
                    clade_key = header
                # You could also warn if clade_field >= len(parts), but
                # at minimum we must return something
            else:
                clade_key = header

            alt_seq = str(record.seq)
            muts: List[str] = []
            for idx, (r_base, a_base) in enumerate(zip(ref_seq, alt_seq)):
                if r_base == "-" or a_base == "-":
                    continue
                if r_base != a_base:
                    pos = ref_map[idx]  # guaranteed non-None here
                    muts.append(f"{r_base}{pos}{a_base}")

            yield clade_key, muts
