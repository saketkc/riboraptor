import numpy as np
import pandas as pd
from .helpers import Interval
from .fasta import FastaReader
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_rna


def counts_to_tpm(counts, sizes):
    """Counts to TPM
    Parameters
    ----------
    counts: array like
            Series/array of counts
    sizes: array like
           Series/array of region sizes
    """
    rate = np.log(counts).subtract(np.log(sizes))
    denom = np.log(np.sum(np.exp(rate)))
    tpm = np.exp(rate - denom + np.log(1e6))
    return tpm


def featurecounts_to_tpm(fc_f, outfile):
    """Convert htseq-counts file to tpm
    Parameters
    ----------
    fc_f: string
             Path to htseq-count output
    outfile: string
             Path to output file with tpm values
    """
    feature_counts = pd.read_csv(fc_f, sep="\t")
    feature_counts = feature_counts.set_index("Geneid")
    feature_counts = feature_counts.drop(columns=["Chr", "Start", "End", "Strand"])
    lengths = feature_counts["Length"]
    feature_counts = feature_counts.drop(columns=["Length"])
    tpm = feature_counts.apply(lambda x: counts_to_tpm(x, lengths), axis=0)
    tpm.columns = [
        col.replace("bams_unique/", "").replace(".bam", "") for col in tpm.columns
    ]
    tpm.to_csv(outfile, sep="\t", index=True, header=True)


def read_ribocounts_raw_count(filepath):
    df = pd.read_csv(filepath, sep="\t")
    df["normalized_count"] = np.divide(1 + df["count"], df["length"])
    df = df.sort_values(by=["gene_id"])
    df = df.set_index("gene_id")
    return df[["count"]]


def read_raw_counts_into_matrix(count_files):
    counts_df = pd.DataFrame()
    for f in count_files:
        filename = path_leaf(f)
        filename = filename.replace("_counts_cnt.txt", "")
        print(filename)
        df = read_ribocounts_raw_count(f)
        df.columns = [filename]
        counts_df = counts_df.join(df, how="outer")
    return counts_df


def ribotricer_index_to_fasta(ribotricer_index, fasta_file):
    """Convert ribotricer index to fasta.

    Parameters
    ----------
    ribotricer_index: str or pd.DataFrame
                      Path to ribotricer index file or read_csv version
    fasta_file: str
                Location of fasta

    Returns
    -------
    sequences: list
               list of list with orf_id and corresponding sequence

    """
    fasta = FastaReader(fasta_file)
    if isinstance(ribotricer_index, str):
        ribotricer_index = pd.read_csv(ribotricer_index, sep="\t").set_index("ORF_ID")

    sequences = []
    for idx, row in ribotricer_index.iterrows():
        intervals = []
        for coordinate in row.coordinate.split(","):
            start, stop = coordinate.split("-")
            start = int(start)
            stop = int(stop)
            interval = Interval(row.chrom, start, stop)
            intervals.append(interval)
        sequence = fasta.query(intervals)
        sequence = "".join(sequence)
        if row.strand == "-":
            sequence = fasta.reverse_complement(sequence)
        sequences.append([idx, sequence])
    return sequences


def get_translated_sequences(sequences):
    for idx, row in enumerate(sequences):
        orf_id, dna_seq = row
        translated_seq = str(Seq(dna_seq, generic_dna).translate())
        sequences[idx].append(translated_seq)
    return sequences


def ribotricer_index_row_to_fasta(ribotricer_index, fasta_file):
    """Convert ribotricer index to fasta (just one row).
    
    Parameters
    ----------
    ribotricer_index: str or pd.DataFrame
                      Path to ribotricer index file or read_csv version
    fasta_file: str
                Location of fasta
    
    Returns
    -------
    sequences: list
               list of list with orf_id and corresponding sequence
    
    """
    fasta = FastaReader(fasta_file)
    # sequences = []
    # idx = row.index

    intervals = []
    for coordinate in row.coordinate.split(","):
        start, stop = coordinate.split("-")
        start = int(start)
        stop = int(stop)
        interval = Interval(row.chrom, start, stop)
        intervals.append(interval)
    sequence = fasta.query(intervals)
    sequence = "".join(sequence)
    if row.strand == "-":
        sequence = fasta.reverse_complement(sequence)
    return sequence
