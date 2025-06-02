### Generating an MSA for Constellation Creation

When using the `alvoc constellations msa` subcommand, you need to supply a multiple‐sequence alignment (MSA) in FASTA format that contains all reference genomes (or representative sequences) from the lineages you wish to capture. Below is a step‐by‐step guide for building such an MSA using common tools (e.g., MAFFT), along with recommendations on sequence selection and basic quality checks.

---

#### 1. Collect Reference Sequences

1. **Decide which genomes to include.**
   Typically, you’ll want:

   - A “reference” genome (e.g., the earliest complete genome for your virus of interest, such as Wuhan‐Hu‐1 for SARS‐CoV‐2).
   - One or more representative sequences from each lineage or clade you care about. These can be downloaded from GenBank or GISAID.
   - Optional: Include outgroup sequences (closely‐related species/strains) if you want to root the tree.

2. **Assemble them into a single FASTA file.**
   For instance, create a file named `references.fasta` containing entries like:

   ```fasta
   >Wuhan-Hu-1|NC_045512.2
   ATTAAAGGTTT... (complete reference sequence)
   >LineageA.23.1_rep1|MW123456.1
   ATGTTGGTTT... (sequence from a representative of A.23.1)
   >LineageB.1.1.7_rep1|MW654321.1
   ATTTAGTGTT... (sequence from B.1.1.7, etc.)
   ```

   Make sure each header is unique and includes enough metadata (accession or lineage name) to identify the sequence later.

---

#### 2. Run MAFFT to Create the Alignment

MAFFT is a widely‐used aligner for viral genomes because it handles large numbers of sequences quickly and accurately.

1. **Install MAFFT** (if not already available). For example, on Ubuntu:

   ```bash
   sudo apt-get update && sudo apt-get install mafft
   ```

   Or via conda:

   ```bash
   conda install -c bioconda mafft
   ```

2. **Execute MAFFT with sensible parameters.**
   A common invocation for full‐length viral genomes is:

   ```bash
   mafft --auto references.fasta > references_aligned.fasta
   ```

   - `--auto` lets MAFFT pick the best algorithm based on your dataset size.
   - If you have many (e.g., >100) long genomes, you might replace `--auto` with `--thread 8 --retree 2 --maxiterate 1000`. For example:

     ```bash
     mafft --thread 8 --retree 2 --maxiterate 1000 references.fasta > references_aligned.fasta
     ```

3. **Inspect the output.**
   Open `references_aligned.fasta` in any alignment viewer (e.g., AliView, Jalview, or even a text editor) and check that:

   - All sequences have roughly the same length.
   - There are no large blocks of gaps at one end (indicative of truncated sequences).
   - Highly divergent sequences still align sensibly (no spuriously long insertions).

---

#### 3. (Optional) Trim or Mask Low‐Confidence Regions

Depending on your downstream analysis and how varied your input set is, you may want to remove poorly aligned regions:

1. **Identify regions with excessive gaps.**
   For example, if more than 50% of sequences have a gap at a given position, you might mask that column.

2. **Use a simple script or tool (e.g., trimAl) to perform automated trimming.**
   Example with trimAl:

   ```bash
   trimal -in references_aligned.fasta -out references_aligned_trimmed.fasta \
         -automated1
   ```

   - `-automated1` applies a default trimming strategy suitable for most alignments.

3. **Double‐check your trimmed FASTA.**
   Ensure that no essential genomic regions (e.g., S gene for SARS‐CoV‐2) were accidentally clipped.

---

#### 4. Validate the Final MSA

Before feeding the MSA into `alvoc constellations msa`, run a quick sanity check:

- **Sequence lengths should be consistent.**

  ```bash
  awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next;} {seqlen += length($0);} END{print seqlen}' \
      references_aligned_trimmed.fasta | sort -nu
  ```

  - All printed lengths should be the same (or very close if you intentionally kept variable‐length UTRs/gaps).

- **Look for unexpected characters.**
  Make sure there are only `A, T, G, C, N, -` (and possibly lowercase) in your FASTA. No spurious symbols like `*` or spaces.

---

#### 5. Run `alvoc constellations msa`

Once you have a validated MSA (e.g., `references_aligned_trimmed.fasta`), call the constellations subcommand. Ensure you specify your header delimeter and clade field index so the algorithm correctly groups your constellations:

```bash
alvoc constellations msa \
    references_aligned_trimmed.fasta \
    -o constellations
    -cd '|' \
    -cf 0
```

After this runs, you’ll get a directory as specified by the output parameter, containing a `constellations.json` file that specifies lineage‐specific site mutations derived from your MSA.

```sh
constellations/
├── constellations.json
├── constellations.json.sha256
├── constellations.manifest.json
└── lineages.txt
```

This JSON can then be passed into `alvoc find-lineages` or `alvoc find-mutants`.

---

#### 6. Tips & Best Practices

- **Use a consistent reference sequence.**
  All coordinates in downstream mutation files are referenced against the same genome. If you align to Wuhan‐Hu-1 (NC_045512.2), make sure all coordinates reported by Alvoc are interpretable relative to that locus.

- **Keep your input FASTA headers simple yet unique.**
  Avoid spaces or special characters; use underscores (`_`) or pipes (`|`) instead. The `--metadata` CSV will match on the exact header line (everything after the `>`), so they must correspond exactly.

- **Document your MSA parameters.**
  If you plan to compare constellations across different runs, note which MAFFT flags (e.g., `--auto` vs.\ `--globalpair`) you used. Small differences in alignment can lead to slight shifts in how site mutations are called.

- **Re‐use alignments whenever possible.**
  If you expect to add new lineages or update an existing constellation, simply append new sequences to your old `references.fasta` and re‐align. You can keep your previous `--metadata` CSV and only add rows for the new headers.

---
