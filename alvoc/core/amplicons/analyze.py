from pathlib import Path
from typing import Callable
import pysam
from alvoc.core.amplicons.visualize import plot_depths, plot_depths_gc
import pandas as pd


class Amplicons:
    """
    Analyzes sample BAM files for coverage and GC content depth correlation.

    Attributes:
        file_path: Path to the file listing samples or a single BAM file.
        inserts: Path to the CSV file detailing the regions (amplicons) to evaluate.
        sequence: DNA sequence used for reference in analysis.
        outdir: Output directory where results will be saved.

    Methods:
        amplicon_coverage(): Determines and plots amplicon coverage.
        gc_depth(): Determines and plots the GC depth correlation.

    Examples:
        Example of a text file listing sample BAM files:
        ```
        /path/to/sample1.bam Sample 1
        /path/to/sample2.bam Sample 2
        ```

        Example of a CSV file detailing regions:
        ```
        chromosome,start,end
        chr1,10000,10500
        chr2,20000,20500
        ```
    """

    def __init__(self, file_path: Path, inserts: Path, sequence: str, outdir: Path):
        self.file_path = file_path
        self.inserts = inserts
        self.sequence = sequence
        self.outdir = outdir
        try:
            self.bed_df = pd.read_csv(inserts)
            if not {"chrom", "chromStart", "chromEnd"}.issubset(self.bed_df.columns):
                raise ValueError(
                    "CSV file must contain 'chrom', 'chromStart', 'chromEnd' columns."
                )
        except pd.errors.EmptyDataError:
            raise FileNotFoundError(
                "Failed to read or validate inserts file: The file is empty or invalid"
            )
        except FileNotFoundError:
            raise FileNotFoundError(
                "Failed to read or validate inserts file: The file was not found"
            )

    @staticmethod
    def calculate_gc_content(sequence: str) -> float:
        """
        Calculate the GC content of a given sequence.

        Args:
            sequence: DNA sequence.

        Returns:
            GC content as a fraction.
        """
        gc_count = sequence.count("G") + sequence.count("C")
        return gc_count / len(sequence) if len(sequence) > 0 else 0

    def amplicon_coverage(self):
        """Determines and plots amplicon coverage for the specified samples."""
        self.process_samples(plot_depths)

    def gc_depth(self):
        """Determines and plots the GC depth correlation for the specified samples."""
        self.process_samples(plot_depths_gc)

    def process_samples(self, plot_function: Callable):
        """
        Processes sample files to compute and plot depth information.

        Args:
            plot_function (Callable): Function used to plot the computed data.
        """
        final_inserts = self._prepare_inserts()
        sample_results, sample_names = [], []
        if self.file_path.suffix == ".bam":
            sample_results.append(
                self.find_depths_in_bam(self.file_path, final_inserts)
            )
            sample_names.append("")
        else:
            with self.file_path.open("r") as f:
                samples = [
                    line.strip().split("\t") for line in f.readlines() if line.strip()
                ]
            for bam, sample in samples:
                bam_path = Path(bam)
                if bam_path.suffix == ".bam":
                    sample_results.append(
                        self.find_depths_in_bam(bam_path, final_inserts)
                    )
                    sample_names.append(sample)
        plot_function(sample_results, sample_names, self.inserts, self.outdir)

    def _prepare_inserts(self):
        """
        Prepares amplicon inserts data with calculated GC content for each region.
        """
        gc_contents = []
        for _, row in self.bed_df.iterrows():
            section = self.sequence[int(row["chromStart"]) : int(row["chromEnd"])]
            gc_contents.append(self.calculate_gc_content(section))
        self.bed_df["gcContent"] = gc_contents
        return self.bed_df.values.tolist()

    def find_depths_in_bam(
        self, bam_path: Path, inserts: list[list], max_depth: int = 50000
    ) -> dict:
        """
        Reads a BAM file and computes the depth of reads at positions defined by `inserts`.

        Args:
            bam_path: Path to the BAM file.
            inserts: List of Lists containing information about the regions (amplicons) of interest.
            max_depth: Maximum depth to be considered to prevent memory overflow.

        Returns:
            A dictionary mapping amplicon identifiers to their corresponding read depth.
        """
        amplified = {}
        with pysam.AlignmentFile(bam_path.as_posix(), "rb") as samfile:
            amp_mids = {int((int(i[1]) + int(i[2])) / 2): i[3] for i in inserts}
            amplified = {i[3]: 0 for i in inserts}
            for pileupcolumn in samfile.pileup(max_depth=max_depth):
                pos = pileupcolumn.reference_pos
                if pos in amp_mids:
                    depth = pileupcolumn.get_num_aligned()
                    amplified[amp_mids[pos]] = depth

        return amplified
