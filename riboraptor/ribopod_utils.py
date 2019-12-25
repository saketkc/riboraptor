import pandas as pd
import numpy as np


def get_mean_median_std_ribotricer(translating_orf_file):
    """Get mean, median, and standard deviation of protein-coding ORFs"""
    df = pd.read_csv(translating_orf_file, sep='\t', usecols= ['phase_score', 'transcript_type'])
    df_subset = df[df.transcript_type == 'protein_coding']
    mean = np.mean(df_subset.phase_score)
    median = np.median(df_subset.phase_score)
    std = np.std(df_subset.phase_score)
    return pd.Series([mean, median, std])

def row_mean_median_std(row):
    """Get mean, median, and standard deviation for a row of translating_ORFs"""
    srp = row.study_accession
    srx = row.experiment_accession
    path = row.srp_path
    filepath = os.path.join(path, 'ribotricer_results', '{}_translating_ORFs.tsv'.format(srx))
    if os.path.exists(filepath):
        return get_mean_median_std_ribotricer(filepath)
    return pd.Series([np.nan, np.nan, np.nan])
