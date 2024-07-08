import pysam
from alvoc.core.amplicons.visualize import plot_depths, plot_depths_gc

def amplicon_coverage(file_path, inserts):
    """
    Determines and plots amplicon coverage for samples listed in a file or a single BAM file.

    Args:
    file_path (str): Path to the file listing samples or a single BAM file.
    inserts (list): List of tuples detailing the regions (amplicons) to evaluate.
    """
    process_samples(file_path, plot_depths, inserts)

def gc_depth(file_path, inserts):
    """
    Determines and plots the GC depth correlation for samples listed in a file or a single BAM file.

    Args:
    file_path (str): Path to the file listing samples or a single BAM file.
    inserts (list): List of tuples detailing the regions (amplicons) to evaluate.
    """
    process_samples(file_path, plot_depths_gc, inserts)

def process_samples(file_path, plot_function, inserts):
    """
    Processes a file containing sample paths and names, extracts depths from BAM files, and visualizes them.

    Args:
    file_path (str): Path to a file listing sample BAM files and their labels or a single BAM file.
    plot_function (function): Function used for plotting results.
    inserts (list): List of tuples detailing the regions (amplicons) to evaluate.
    """
    sample_results, sample_names = [], []
    if file_path.endswith('.bam'):
        sample_results.append(find_depths_in_bam(file_path, inserts))
        sample_names.append('')
    else:
        with open(file_path, 'r') as f:
            samples = [line.strip().split('\t') for line in f.readlines() if line.strip()]
        for sample in samples:
            if sample[0].endswith('.bam'):
                sample_results.append(find_depths_in_bam(sample[0], inserts))
                sample_names.append(sample[1])
    plot_function(sample_results, sample_names, inserts)

def find_depths_in_bam(bam_path, inserts, max_depth=50000):
    """
    Reads a BAM file and computes the depth of reads at positions defined by `inserts`.

    Args:
    bam_path (str): Path to the BAM file.
    inserts (list): List of tuples containing information about the regions (amplicons) of interest.
    max_depth (int): Maximum depth to be considered to prevent memory overflow.

    Returns:
    dict: A dictionary mapping amplicon identifiers to their corresponding read depth.
    """
    samfile = pysam.AlignmentFile(bam_path, "rb")
    amp_mids = {int((int(i[1]) + int(i[2])) / 2): i[3] for i in inserts}
    amplified = {i[3]: 0 for i in inserts}
    for pileupcolumn in samfile.pileup(max_depth=max_depth):
        pos = pileupcolumn.reference_pos
        if pos in amp_mids:
            depth = pileupcolumn.get_num_aligned() 
            amplified[amp_mids[pos]] = depth

    samfile.close()
    return amplified
