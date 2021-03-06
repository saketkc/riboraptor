import pysam
from .helpers import is_read_uniq_mapping, create_bam_index
from .helpers import read_refseq_bed, read_bed_as_intervaltree
from collections import Counter
import numpy as np


def infer_protocol(bam, bed, n_reads=500000, drop_probability=0.2):
    """Infer strandedness protocol given a bam file

    Parameters
    ----------
    bam: string
         Path to bam file
    bed: string
         Path to gene bed file
    n_reads: int
             Number of reads to use (downsampled)
    drop_probability: float [0-1]
                      Reads are randomly droppped with this probability to
                      simulate sampling


    Returns
    -------
    protocol: string
              unstranded/forward/reverse
    forward_mapped_reads: float
          Proportion of reads of type + mapping to + (++) or - mapping to - (--)
    reverse_mapped_reads: float
          Proportion of reads of type + mapping to - (+-) or - mapping to + (-+)
    """
    np.random.seed(42)
    iteration = 0
    create_bam_index(bam)
    bam = pysam.AlignmentFile(bam, "rb")
    bed = read_bed_as_intervaltree(bed)
    strandedness = Counter()
    for read in bam.fetch():
        if drop_probability > 0:
            should_drop = np.random.choice(
                [1, 0], p=[drop_probability, 1 - drop_probability]
            )
            if should_drop:
                continue
        if not is_read_uniq_mapping(read):
            continue
        if read.is_reverse:
            mapped_strand = "-"
        else:
            mapped_strand = "+"
        mapped_start = read.reference_start
        mapped_end = read.reference_end
        chrom = read.reference_name
        gene_strand = list(set(bed[chrom].find(mapped_start, mapped_end)))
        if len(gene_strand) != 1:
            # Filter out genes with ambiguous strand info
            # (those) that have a tx_start on opposite strands
            continue
        gene_strand = gene_strand[0]
        strandedness["{}{}".format(mapped_strand, gene_strand)] += 1
        iteration += 1
        if iteration >= n_reads:
            break
    ## Add pseudocounts
    strandedness["++"] += 1
    strandedness["--"] += 1
    strandedness["+-"] += 1
    strandedness["-+"] += 1

    total = sum(strandedness.values())
    forward_mapped_reads = (strandedness["++"] + strandedness["--"]) / total
    reverse_mapped_reads = (strandedness["-+"] + strandedness["+-"]) / total
    ratio = forward_mapped_reads / reverse_mapped_reads
    # Prefer checking for unstrandedness
    # Check if the forward mapped reads - 0.5 is small,
    # this threhold is defined to be 0.1
    if np.isclose([np.abs(forward_mapped_reads - 0.5)], [0], atol=0.1):
        return "unstranded", forward_mapped_reads, reverse_mapped_reads, total
    elif forward_mapped_reads >= 0.5:
        return "forward", forward_mapped_reads, reverse_mapped_reads, total
    else:
        return "reverse", forward_mapped_reads, reverse_mapped_reads, total
