import numpy as np
import pandas as pd


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
