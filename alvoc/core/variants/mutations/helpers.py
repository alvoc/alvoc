from pysam import PileupColumn


def mut_in_col(pileupcolumn: PileupColumn, mut: str):
    """Count the occurrences and non-occurrences of a mutation in a pileup column.

    Args:
        pileupcolumn : A pileup column from a BAM file.
        mut : The mutation to count.

    Returns:
        tuple: A tuple containing counts of mutations and non-mutations.
    """
    muts = not_muts = 0
    for pileupread in pileupcolumn.pileups:
        qpos = pileupread.query_position
        if qpos is None:
            not_muts += 1
            continue
        base = pileupread.alignment.query_sequence
        if base is not None:
            base = base[qpos]
        if base == mut:
            muts += 1
        else:
            not_muts += 1
    return muts, not_muts


def print_mut_results(mut_results, min_depth):
    """Print the results of mutation detection.

    Args:
        mut_results (dict): Dictionary containing mutation results.
        min_depth (int): Minimum depth to consider for reporting.
    """
    cov = mut_cov = 0
    for name, (muts, not_muts) in mut_results.items():
        total = muts + not_muts
        if total >= min_depth:
            cov += 1
            if muts > 0:
                mut_cov += 1

    print("{}/{} mutations covered".format(cov, len(mut_results)))
    print("{}/{} mutations detected".format(mut_cov, len(mut_results)))
