import os
import pandas as pd
import numpy as np


def get_mean_median_std_ribotricer(translating_orf_file, orf_type):
    """Get mean, median, and standard deviation of protein-coding ORFs"""
    if isinstance(orf_type, str):
        orf_type = [orf_type]
    df = pd.read_csv(
        translating_orf_file,
        sep="\t",
        usecols=["phase_score", "ORF_type", "transcript_type"],
    )
    # df_subset = df[df.transcript_type == 'protein_coding']
    df_subset = df[df.ORF_type.isin(orf_type)]
    mean = np.mean(df_subset.phase_score)
    median = np.median(df_subset.phase_score)
    std = np.std(df_subset.phase_score)
    return pd.Series([mean, median, std])


def row_mean_median_std(row, orf_type):
    """Get mean, median, and standard deviation for a row of translating_ORFs"""
    srp = row.study_accession
    srx = row.experiment_accession
    path = row.srp_path
    filepath = os.path.join(
        path, "ribotricer_results", "{}_translating_ORFs.tsv".format(srx)
    )
    if os.path.exists(filepath):
        return get_mean_median_std_ribotricer(filepath, orf_type)
    return pd.Series([np.nan, np.nan, np.nan])
