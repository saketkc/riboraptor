from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import tqdm
import numpy as np
import re
import pickle
import os

import pandas as pd
from .helpers import mkdir_p
from .helpers import symlink_force


def load_tpm(path):
    df = pd.read_table(path, names=['gene_id', 'tpm']).set_index('gene_id')
    return df


def get_cell_line_or_tissue(row):
    if str(row['cell_line']).strip() and str(
            row['cell_line']).strip() != 'nan':
        return '{}-{}-{}'.format(row['cell_line'], row['study_accession'],
                                 row['experiment_accession'])
    if str(row['tissue']).strip() and str(row['tissue']).strip() != 'nan':
        return '{}-{}-{}'.format(row['tissue'], row['study_accession'],
                                 row['experiment_accession'])
    if str(row['source_name']).strip() and str(
            row['source_name']).strip() != 'nan':
        return '{}-{}-{}'.format(row['source_name'], row['study_accession'],
                                 row['experiment_accession'])
    if row['study_accession'].strip() == 'SRP052229':
        print(row)
    return '{}-{}-{}'.format(row['source_name'], row['study_accession'],
                             row['experiment_accession'])


def determine_cell_type(sample_attribute):
    sample_attribute = str(sample_attribute)
    if 'cell line:' in sample_attribute:
        x = re.search(r'cell line: \w+', sample_attribute)
        return x.group(0).strip('cell line: ').rstrip(' ').upper()
    if 'cell_line:' in sample_attribute:
        x = re.search(r'cell_line: \w+', sample_attribute)
        return x.group(0).strip('cell_line: ').rstrip(' ').upper()
    if 'cell-line:' in sample_attribute:
        x = re.search(r'cell-line: \w+', sample_attribute)
        return x.group(0).strip('cell-line: ').rstrip(' ').upper()
    if 'cell_type:' in sample_attribute:
        x = re.search(r'cell_type: \w+', sample_attribute)
        return x.group(0).strip('cell_type: ').rstrip(' ').upper()
    if 'source_name:' in sample_attribute:
        x = re.search(r'source_name: \w+', sample_attribute)
        return x.group(0).strip('source_name: ').rstrip(' ').upper()
    else:
        #pass
        print('Found {}'.format(sample_attribute))
        return np.nan


def get_tissue_type(sample_attribute):
    sample_attribute = str(sample_attribute)
    if 'tissue: ' in sample_attribute:
        x = re.search(r'tissue: \w+', sample_attribute)
        return x.group(0).strip('tissue: ').rstrip(' ').lower()
    else:
        print('Found {}'.format(sample_attribute))
        return np.nan


def get_strain_type(sample_attribute):
    sample_attribute = str(sample_attribute)
    if 'strain: ' in sample_attribute:
        x = re.search(r'strain: \w+', sample_attribute)
        return x.group(0).strip('strain: ').rstrip(' ').lower()
    else:
        print('Found {}'.format(sample_attribute))
        return np.nan


def summary_starlogs_over_runs(directory, list_of_srr):
    df = pd.DataFrame()
    files_not_found = []
    for run in list_of_srr:
        if not os.path.isfile(os.path.join(directory, run + '.tsv')):
            files_not_found.append(run)
            continue
        temp_df = pd.read_table(os.path.join(directory, run + '.tsv'))
        df = pd.concat([df, temp_df])
    return df, files_not_found


def get_enrichment_cds_stats(pickle_file):
    data = pickle.load(open(pickle_file, 'rb'))
    mean = np.nanmean(data.values())
    median = np.nanmedian(data.values())
    stddev = np.nanstd(data.values())
    minx = np.nanmin(data.values())
    maxx = np.nanmax(data.values())
    return minx, maxx, mean, median, stddev


def get_fragment_enrichment_score(txt_file):
    with open(txt_file) as fh:
        data = fh.read()
    enrichment = data.strip('\(').strip('\)').strip(' ').strip()
    enrichment, pval = enrichment.split(',')
    if 'nan' not in enrichment:
        return float(enrichment.strip('Enrichment: ')), float(
            pval.strip(')').strip('pval: '))
    else:
        return np.nan, 1


def filter_single_end_samples(df):
    """Filter single end samples from a dataframe

    Parameters
    ----------
    df: DataFrame
        Dataframe as obtained from SRAb.sra_convert()

    Returns
    -------
    df: DataFrame
        DataFrame with only single end samples
    """
    df = df[~df['library_strategy'].str.contains('PAIRED')]
    return df


def copy_sra_data(
        df,
        taxon_id_map={
            10090: 'mm10',
            9606: 'hg38'
        },
        sra_source_dir='/staging/as/skchoudh/SRA_datasets/',
        sra_dest_dir='/staging/as/skchoudh/re-ribo-datasets/samples_to_process/'
):
    """Copy SRA data to a new location retaining only single ended samples."""
    df = filter_single_end_samples(df)
    assert len(df.study_accession.unique()) == 1, 'Multiple SRPs found'
    srp = df.study_accession.unique()[0]
    df_grouped = df.groupby(['taxon_id'])
    srp_source_dir = os.path.join(sra_source_dir, srp)

    for taxon_id, df_group in df_grouped:
        species = taxon_id_map[taxon_id]
        species_dest_dir = os.path.join(sra_dest_dir, species)
        srp_dest_dir = os.path.join(species_dest_dir, srp)
        mkdir_p(os.path.join(species_dest_dir, srp))
        source_loc = srp_source_dir + os.path.sep + df[
            'experiment_accession'].str.cat(
                df['run_accession'] + '.sra', sep=os.path.sep)
        dest_loc = srp_dest_dir + os.path.sep + df[
            'experiment_accession'].str.cat(
                df['run_accession'] + '.sra', sep=os.path.sep)
        for source, dest in tqdm(zip(source_loc, dest_loc)):
            mkdir_p(os.path.dirname(dest))
            if os.path.isfile(source):
                symlink_force(source, dest)
