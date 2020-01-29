import os
import pandas as pd
import numpy as np
from .helpers import mkdir_p
from .helpers import path_leaf
from .ribotricer_utils import read_raw_counts_into_matrix


def get_mean_median_std_ribotricer_orftype(translating_orf_file, orf_type):
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


def row_mean_median_std_ribotricer_orftype(row, orf_type):
    """Get mean, median, and standard deviation for a row of translating_ORFs"""
    srx = row.experiment_accession
    path = row.srp_path
    filepath = os.path.join(
        path, "ribotricer_results", "{}_translating_ORFs.tsv".format(srx)
    )
    if os.path.exists(filepath):
        return get_mean_median_std_ribotricer_orftype(filepath, orf_type)
    return pd.Series([np.nan, np.nan, np.nan])


def get_mean_median_std_ribotricer_txtype(translating_orf_file, txtype):
    """Get mean, median, and standard deviation of given transcript type

    txtype: list

    """
    if isinstance(txtype, str):
        txtype = [txtype]
    df = pd.read_csv(
        translating_orf_file,
        sep="\t",
        usecols=["phase_score", "ORF_type", "transcript_type"],
    )
    df_subset = df[df.transcript_type.isin(txtype)]
    mean = np.mean(df_subset.phase_score)
    median = np.median(df_subset.phase_score)
    std = np.std(df_subset.phase_score)
    return pd.Series([mean, median, std])


def get_gene_mean_median_max_std_ribotricer(translating_orf_file, only_ATG=True):
    """Get mean, median, and standard deviation of given transcript type

    txtype: list

    """
    df = pd.read_csv(
        translating_orf_file,
        sep="\t",
        usecols=["phase_score", "ORF_type", "start_codon", "gene_id"],
    )
    if only_ATG:
        df = df.loc[df.start_codon == "ATG"]
    df_subset = df.loc[df.ORF_type == "annotated"]
    df_subset = df_subset.groupby("gene_id").agg(
        min_phase_score=pd.NamedAgg(column="phase_score", aggfunc="min"),
        max_phase_score=pd.NamedAgg(column="phase_score", aggfunc="max"),
        median_phase_score=pd.NamedAgg(column="phase_score", aggfunc="median"),
        mean_phase_score=pd.NamedAgg(column="phase_score", aggfunc="mean"),
        std_phase_score=pd.NamedAgg(column="phase_score", aggfunc="std"),
        n_orfs=pd.NamedAgg(column="phase_score", aggfunc="count"),
    )
    return df_subset


def row_mean_median_std_ribotricer_txtype(row, txtype):
    """Get mean, median, and standard deviation for a row of translating_ORFs"""
    srx = row.experiment_accession
    path = row.srp_path
    filepath = os.path.join(
        path, "ribotricer_results", "{}_translating_ORFs.tsv".format(srx)
    )
    if os.path.exists(filepath):
        return get_mean_median_std_ribotricer_txtype(filepath, txtype)
    return pd.Series([np.nan, np.nan, np.nan])


def collate_ribocounts_in_df(srp_df, region_type):
    """Collect all ribocount files inside a directory and write them as a dataframe"""
    srp = srp_df.study_accession.tolist()[0]
    srp_assembly_grouped = srp_df.groupby("assembly")
    for assembly, srp_assembly_df in srp_assembly_grouped:
        srp_path = srp_assembly_df.srp_path.tolist()[0]
        experiments = list(sorted(srp_assembly_df.experiment_accession.unique()))
        files = [
            os.path.join(
                srp_path,
                "ribotricer_results",
                "region_counts",
                region_type,
                "{}_counts_cnt.txt".format(experiment),
            )
            for experiment in experiments
        ]
        files = [f for f in files if os.path.exists(f)]
        if not files:
            continue
        collected_counts = read_raw_counts_into_matrix(files)
        collected_counts.columns = [
            "{}_{}_{}_{}".format(assembly, srp, srx, region_type)
            for srx in collected_counts.columns
        ]
        collected_counts = collected_counts.reset_index()
        dir_to_write = os.path.join(
            srp_path, "ribotricer_results", "ribotricer_ribo_counts_df"
        )
        mkdir_p(dir_to_write)
        file_to_write = os.path.join(
            dir_to_write, "{}_{}_{}.tsv".format(assembly, srp, region_type)
        )
        collected_counts.to_csv(file_to_write, sep="\t", index=False, header=True)


def aggregate_ribotricer_output_to_gene(infile, outfile):
    """Aggregate counts so that they are represented as one record per gene"""
    df = pd.read_csv(
        infile,
        sep="\t",
        usecols=[
            "ORF_type",
            "gene_id",
            "length",
            "valid_codons",
            "read_count",
            "phase_score",
        ],
    )
    df["length_normalized_read_count"] = df["read_count"].divide(df["length"])
    df["length_normalized_read_count"] = np.round(df["length_normalized_read_count"], 3)
    df["valid_codon_ratio"] = df["valid_codons"].divide(df["length"])
    df["valid_codon_ratio"] = np.round(df["valid_codon_ratio"], 3)
    df["phase_score"] = np.round(df["phase_score"], 3)
    df = df[
        [
            "gene_id",
            "ORF_type",
            "valid_codon_ratio",
            "length_normalized_read_count",
            "phase_score",
        ]
    ]
    df_grouped = df.groupby(["gene_id", "ORF_type"]).agg(list)
    df_grouped = df_grouped.reset_index()
    mkdir_p(os.path.dirname(outfile))

    for region_type in df_grouped.ORF_type.unique():
        df_subset = df_grouped.loc[df_grouped.ORF_type == region_type]
        df_subset.to_csv(
            os.path.splitext(outfile)[0] + "__{}.tsv".format(region_type),
            sep="\t",
            header=True,
            index=False,
        )
    df_grouped.to_csv(outfile, sep="\t", header=True, index=False)


def read_ribotricer_output(filepath):
    experiment = path_leaf(filepath)
    experiment = experiment.replace("_translating_ORFs.tsv", "")
    df = pd.read_csv(
        filepath,
        sep="\t",
        usecols=["ORF_ID", "phase_score", "valid_codons", "length", "read_count"],
    ).set_index("ORF_ID")
    df["valid_codon_ratio"] = df["valid_codons"].divide(df["length"])

    df_phase_score = df[["phase_score"]].rename(columns={"phase_score": experiment})
    df_vcr = df[["valid_codon_ratio"]].rename(columns={"valid_codon_ratio": experiment})
    df_rc = df[["read_count"]].rename(columns={"read_count": experiment})
    return df_phase_score, df_vcr, df_rc


def read_phase_score_vcr_rc(files):
    phase_score_df = pd.DataFrame()
    vcr_df = pd.DataFrame()
    rc_df = pd.DataFrame()
    for filepath in files:
        p_df, v_df, r_df = read_ribotricer_output(filepath)
        phase_score_df = phase_score_df.join(p_df)
        vcr_df = vcr_df.join(v_df)
        rc_df = rc_df.join(r_df)
    return phase_score_df, vcr_df, rc_df


def _rename_cols(collected_df, assembly, srp):
    collected_df.columns = [
        "{}_{}_{}".format(assembly, srp, srx) for srx in collected_df.columns
    ]
    collected_df = collected_df.reset_index()
    return collected_df


def collect_ribotricer_orfs(srp_df):
    srp = srp_df.study_accession.tolist()[0]
    srp_assembly_grouped = srp_df.groupby("assembly")
    for assembly, srp_assembly_df in srp_assembly_grouped:
        srp_path = srp_assembly_df.srp_path.tolist()[0]
        experiments = list(sorted(srp_assembly_df.experiment_accession.unique()))
        files = [
            os.path.join(
                srp_path,
                "ribotricer_results",
                "{}_translating_ORFs.tsv".format(experiment),
            )
            for experiment in experiments
        ]

        files = [f for f in files if os.path.exists(f)]
        if not files:
            continue

        phase_score_df, vcr_df, rc_df = read_phase_score_vcr_rc(files)

        dir_to_write = os.path.join(
            srp_path, "ribotricer_results", "ribotricer_orfs_df"
        )
        mkdir_p(dir_to_write)
        file_to_write = os.path.join(dir_to_write, "{}_{}".format(assembly, srp))
        phase_score_df = _rename_cols(phase_score_df, assembly, srp)
        vcr_df = _rename_cols(vcr_df, assembly, srp)
        rc_df = _rename_cols(rc_df, assembly, srp)

        phase_score_df.to_csv(
            file_to_write + "_phase_score.tsv", sep="\t", index=False, header=True
        )
        vcr_df.to_csv(
            file_to_write + "_phase_score.tsv", sep="\t", index=False, header=True
        )
        rc_df.to_csv(
            file_to_write + "_read_count.tsv", sep="\t", index=False, header=True
        )
