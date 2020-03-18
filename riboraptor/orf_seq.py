import numpy as np
import pyfaidx
import sys
from tqdm import tqdm

import pandas as pd
from .interval import Interval
from .fasta import FastaReader


def offset_start_stop(
    start, stop, strand, chrom_size, upstream_5p_offset, downstream_3p_offset
):
    if strand == "+":
        start = max(start - upstream_5p_offset, 1)
        stop = min(stop + downstream_3p_offset, chrom_size)
    elif strand == "-":
        start = max(start - downstream_3p_offset, 1)
        stop = min(stop + upstream_5p_offset, chrom_size)
    return start, stop


def orf_seq(
    ribotricer_index,
    genome_fasta,
    saveto,
    upstream_5p_offset=0,
    downstream_3p_offset=0,
    translate=False,
):
    """Generate sequence for ribotricer annotation.

  Parameters
  -----------
  ribotricer_index: string
                         Path to ribotricer generate annotation
  genome_Fasta: string
                Path to genome fasta
  saveto: string
          Path to output
  """
    fasta = FastaReader(genome_fasta)
    annotation_df = pd.read_csv(ribotricer_index, sep="\t")
    with open(saveto, "w") as fh:
        fh.write("ORF_ID\tsequence\n")
        for idx, row in tqdm(annotation_df.iterrows(), total=annotation_df.shape[0]):
            chrom = str(row.chrom)
            orf_id = row.ORF_ID
            coordinates = row.coordinate.split(",")
            strand = row.strand
            try:
                chrom_size = fasta.chromosomes[chrom]
            except:
                chrom_size = np.inf
            intervals = []
            seq = ""
            for index, coordinate in enumerate(coordinates):
                start, stop = coordinate.split("-")
                start = int(start)
                stop = int(stop)
                if index == 0:
                    start, stop = offset_start_stop(
                        start,
                        stop,
                        strand,
                        chrom_size,
                        upstream_5p_offset,
                        downstream_3p_offset,
                    )
                interval = Interval(chrom, start, stop, strand)
                intervals.append(interval)

            seq = ("").join(fasta.query(intervals))
            if strand == "-":
                seq = fasta.reverse_complement(seq)
            if translate:
                if len(seq) % 3 != 0:
                    sys.stderr.write(
                        "WARNING: Sequence length with ORF ID '{}' is not a multiple of three. Output sequence might be truncated.\n".format(
                            orf_id
                        )
                    )
                    seq = seq[0 : (len(seq) // 3) * 3]
                seq = translate_nt_to_aa(seq)
            fh.write("{}\t{}\n".format(orf_id, seq))
