---
hide:
  - navigation
---

# Usage

This guide explains how to use **Alvoc**, highlighting key features, input formats, and examples. It assumes you've already installed the tool. If not, please the [installation](index.md) instructions.

---

With Alvoc, you can quickly move from raw sequencing data to actionable insights about viral mutations and lineages.

## Pre-requisites

Alvoc, like many other downstream analysis tools, requires you to have some upstream data to unlock its functionality. There are 2 main requirements for the core commands:

1. A `GB` (genbank) file for a complete genome of the virus to analyze, or an [Entrez](https://www.ncbi.nlm.nih.gov/books/NBK25501/) api search term. Please see the [Extracting Gene Data section](#extracting-gene-data) for a brief on how this works.

2. A sorted `BAM` file, or a `CSV` samplesheet listing samples and absolute paths to the BAM files. An example csv file:

   | sample  | bam                           |
   | ------- | ----------------------------- |
   | sample1 | /absolute/path/to/sample1.bam |
   | sample2 | /absolute/path/to/sample2.bam |

## First steps

Alvoc comes with a variety of tools to help you with abundance learning. Take a gander by running the `help` flag to see a general help message:

<!-- termynal -->

```console
$ alvoc --help
 Usage: alvoc [OPTIONS] COMMAND [ARGS]...                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                            
 Abundance learning for variants of concern                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                            
╭─ Options ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --version                     Current version of Alvoc                                                                                                                                                                                                   │
│ --install-completion          Install completion for the current shell.                                                                                                                                                                                  │
│ --show-completion             Show completion for the current shell, to copy it or customize the installation.                                                                                                                                           │
│ --help                        Show this message and exit.                                                                                                                                                                                                │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Commands ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ find-lineages       Find lineages in samples                                                                                                                                                                                                             │
│ find-mutants        Find mutations in sequencing data, either from BAM files or a sample list.                                                                                                                                                           │
│ extract-gene-data   Extracts gene coordinates and the genome sequence from a GenBank file or generates them using an Entrez API search term.                                                                                                             │
│ amplicons           Get amplicon metrics such as coverage, gc_content and visualizations                                                                                                                                                                 │
│ convert             Tools to convert mutations                                                                                                                                                                                                           │
│ constellations      Tools to make constellations                                                                                                                                                                                                         │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

```

## Extracting Gene Data

Before we take a look at some of the core commands, it is important to understand some of Alvoc's internals. Each command requires gene coordinates and sequence data. We can extract this from genbank files, or we can auto-fetch the data for you using an [Entrez](https://www.ncbi.nlm.nih.gov/books/NBK25501/) search term.

### A brief on the Entrez API

The Entrez API (E-utilities) is a tool provided by the NCBI to interact with its biological databases, such as GenBank and PubMed. It allows you to search for IDs, retrieve data like sequences, and link related records across databases. We use this internally to download a genbank file programatically.

You can test how it works using the following examples to get a feel for it:

1. Access by accession id:

   <!-- termynal -->

   ```console
   $ alvoc extract-gene-data NC_045512.2
   ```

2. You can also just search by taxid, although you need to pass it like this:

   <!-- termynal -->

   ```console
   $ alvoc extract-gene-data  "txid2697049[Organism:exp]"
   ```

3. Sometimes you want to be even more specific

   <!-- termynal -->

   ```console
   $ alvoc extract-gene-data  "txid2697049[Organism:exp] AND \"Wuhan-Hu-1\" AND complete genome"
   ```

!!! note

    Every command in alvoc requires this input or alternatively a genbank file. In general, we recommend using the accession ID for data retrieval, as it is the most straightforward and unambiguous approach. However, the Entrez API is subject to rate limiting, so if you plan on performing multiple successive analyses, we recommend downloading the desired genbank file and supplying it in place of the accession ID.

---

## Finding Lineages

The **`find-lineages`** command identifies viral lineages from sequencing data. Running it will output a lineage abundance heatmap and two csv files. To get started, you'll need your reference genome, a samplesheet listing your bam files, and a constellations file. See the [constellations section for more details](constellations.md).

### Command Example

<!-- termynal -->

```console
$ alvoc find-lineages some_accession_id samples.csv constellations.json --outdir ./results
```
 

## Finding Mutations

The **`find-mutants`** command identifies specific mutations from sequencing and outputs a csv file containing occurences and non-occurences of each mutation, alongside a heatmap.

### Command Example

<!-- termynal -->

```console
$ alvoc find-mutants virus.gb  samples.csv constellations.json mutations.txt --outdir ./results
```

Here's what the mutations.txt file might look like

```text
S:N501Y
G23012A
```

---

## Amplicon Metrics

The **`amplicons`** command analyzes amplicon metrics like coverage and GC content. You'll need to provide:

A `BED` file, containing your inserts of interest. For example:

```text
MN908947.3	50	408	SARS-CoV-2_INSERT_1	1	+
MN908947.3	344	705	SARS-CoV-2_INSERT_2	2	+
MN908947.3	666	1017	SARS-CoV-2_INSERT_3	1	+
```

### Command Example

<!-- termynal -->

```console
$ alvoc amplicons virus.gb samples.csv inserts.bed --outdir ./results
```

---

# Tools

### Converting Mutations

The **`convert`** command transforms mutations between formats.

#### Amino Acid to Nucleotide

<!-- termynal -->

```console
$ alvoc convert aa virus S:N501Y
```

#### Nucleotide to Amino Acid

<!-- termynal -->

```console
$ alvoc convert nt virus G23012A
```
