import pandas
import pysam
from pathlib import Path
import typer
import matplotlib.pyplot as plt

def qc(samples: Path, out: Path):
    # Build a list of (sample_name, bam_path) tuples
    sample_list = []
    if samples.is_dir():
        for bam_path in samples.glob("*.bam"):
            sample_list.append((bam_path.stem, bam_path))
    elif samples.suffix.lower() in (".csv", ".tsv"):
        df = pandas.read_csv(samples)
        if not {"sample", "bam"}.issubset(df.columns):
            typer.echo("CSV must contain 'sample' and 'bam' columns.", err=True)
            raise typer.Exit(code=1)
        for _, row in df.iterrows():
            sample_list.append((row["sample"], Path(row["bam"])))
    else:
        sample_list.append((samples.stem, samples))

    # Collect alignment stats
    records = []
    for sample_name, bam_path in sample_list:
        try:
            samfile = pysam.AlignmentFile(str(bam_path), "rb")
        except Exception as e:
            typer.echo(f"Failed to open BAM {bam_path}: {e}", err=True)
            continue

        total = mapped = proper = duplicates = 0
        for read in samfile.fetch(until_eof=True):
            total += 1
            if read.is_unmapped:
                continue
            mapped += 1
            if read.is_proper_pair:
                proper += 1
            if read.is_duplicate:
                duplicates += 1
        samfile.close()

        pct_mapped = (mapped / total * 100) if total > 0 else 0.0
        records.append(
            {
                "sample": sample_name,
                "total_reads": total,
                "mapped_reads": mapped,
                "properly_paired": proper,
                "duplicates": duplicates,
                "pct_mapped": round(pct_mapped, 2),
            }
        )

    # Output results
    result_df = pandas.DataFrame.from_records(records)
    if out:
        result_df.to_csv(out, index=False)
        typer.echo(f"QC summary written to {out}")

    # Plot and save the histogram of mapping rates
    plot_file = Path(out).with_suffix(".png")
    plot_mapping_rate_hist(result_df, plot_file)
    typer.echo(f"Mapping rate histogram saved to {plot_file}")

def plot_mapping_rate_hist(df, output_file):
    """
    Plot a bar chart of mapping rate per sample with y-axis scaled to 0-100% and save it.
    """
    plt.figure(figsize=(12, 6))
    plt.bar(df['sample'], df['pct_mapped'], color='orange')
    plt.ylim(0, 100)  # ensure full 0–100% scale
    plt.xlabel('Sample')
    plt.ylabel('Mapping Rate (%)')
    plt.title('Mapping Rate per Sample')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()