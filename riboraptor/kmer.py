from collections import Counter
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from tqdm import tqdm
import gzip
import pandas as pd


def fastq_kmer_histogram(fastq_file,
                         kmer_length=range(5, 31),
                         five_prime=False,
                         max_seq=1000000,
                         offset=0):
    """Get a histogram of kmers from a fastq file

    Parameters
    ----------
    fastq_file: string
                Location of .fastq or .fastq.gz
    kmer_length: range
                 Range of kmer to consider
    five_prime: bool
                Should consider sequences from 5' end?
                Default: False (uses sequence from 3' end)
    max_seq: int
             Maximum number of sequences to consider
    offset: int
            Offset to ignore at the 5' end or 3'end
            Example: If the offset is 3, the first 3 bases will be skipped
            and kmers will start only from the 4th position
            For 3'end if the offset is 3, the last 3 bases will be skipped

    Returns
    -------
    kmers: Series
           Sorted series with most frequent kmer

    """
    cur_count = 0
    should_continue = True
    if '.gz' in fastq_file:
        # Open as a gzip file
        handle = gzip.open(fastq_file, 'rt')
    else:
        handle = open(fastq_file, 'r')
    histogram = {k: Counter() for k in kmer_length}

    with tqdm(total=max_seq) as pbar:
        for title, seq, qual in FastqGeneralIterator(handle):
            if not should_continue:
                break
            cur_count += 1
            for k in kmer_length:
                if not five_prime:
                    if not offset:
                        k_seq = seq[-k:]
                    else:
                        k_seq = seq[-k - offset:-offset]
                else:
                    k_seq = seq[offset:k + offset]
                histogram[k][k_seq] += 1
                if cur_count >= max_seq:
                    should_continue = False
            pbar.update()
    handle.close()
    kmers = {}
    for k, v in histogram.items():
        kmers[k] = pd.Series(v).sort_values(ascending=False) / max_seq * 100
    return kmers
