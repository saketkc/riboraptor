from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from tqdm import tqdm
import numpy as np
import re
import pickle
import os

import pandas as pd
from textwrap import dedent
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


re_ribo_root_dir = '/staging/as/skchoudh/SRA_datasets/'
samples_to_process_dir = '/staging/as/skchoudh/re-ribo-datasets/'
re_ribo_config_dir = '/home/cmb-panasas2/skchoudh/github_projects/re-ribo-smk/snakemake/configs'
re_ribo_analysis_dir = '/staging/as/skchoudh/re-ribo-analysis/'
riboraptor_annotation_dir = '/home/cmb-panasas2/skchoudh/github_projects/riboraptor/riboraptor/annotation/'

genome_annotation_map = {
    'hg38': 'v25',
    'mm10': 'vM11',
    'mg1655': '',
    'sacCerR64': 'v91',
    'MG1655': 'ASM584v2.38',
    'BDGP6': 'v91',
    'GRCz10': 'v91'
}
taxon_id_map = {
    10090: 'mm10',
    9606: 'hg38',
    4932: 'sacCerR64',
    511145: 'MG1655',
    7227: 'BDGP6',
    7955: 'GRCz10'
}

genome_fasta_map = {
    'hg38':
    '/home/cmb-panasas2/skchoudh/genomes/hg38/fasta/hg38.fa',
    'mm10':
    '/home/cmb-panasas2/skchoudh/genomes/mm10/fasta/mm10.fa',
    'sacCerR64':
    '/home/cmb-panasas2/skchoudh/genomes/sacCerR64/fasta/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa',
    'MG1655':
    '/home/cmb-panasas2/skchoudh/genomes/escherichia_coli_str_k_12_substr_mg1655/fasta/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa',
    'BDGP6':
    '/home/cmb-panasas2/skchoudh/genomes/drosophila_melanogaster_BDGP6/fasta/Drosophila_melanogaster.BDGP6.dna.toplevel.fa',
    'GRCz10':
    '/home/cmb-panasas2/skchoudh/genomes/GRCz10/fasta/Danio_rerio.GRCz10.dna.toplevel.fa'
}

chrom_sizes_map = {
    'hg38':
    '/home/cmb-panasas2/skchoudh/genomes/hg38/fasta/hg38.chrom.sizes',
    'mm10':
    '/home/cmb-panasas2/skchoudh/genomes/mm10/fasta/mm10.chrom.sizes',
    'sacCerR64':
    '/home/cmb-panasas2/skchoudh/genomes/sacCerR64/fasta/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.sizes',
    'MG1655':
    '/home/cmb-panasas2/skchoudh/genomes/escherichia_coli_str_k_12_substr_mg1655/fasta/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.sizes',
    'BDGP6':
    '/home/cmb-panasas2/skchoudh/genomes/drosophila_melanogaster_BDGP6/fasta/Drosophila_melanogaster.BDGP6.dna.toplevel.sizes',
    'GRCz10':
    '/home/cmb-panasas2/skchoudh/genomes/GRCz10/fasta/Danio_rerio.GRCz10.dna.toplevel.sizes'
}

star_index_map = {
    'hg38':
    '/home/cmb-panasas2/skchoudh/genomes/hg38/star_annotated',
    'mm10':
    '/home/cmb-panasas2/skchoudh/genomes/mm10/star_annotated',
    'sacCerR64':
    '/home/cmb-panasas2/skchoudh/genomes/sacCerR64/star_annotated',
    'MG1655':
    '/home/cmb-panasas2/skchoudh/genomes/escherichia_coli_str_k_12_substr_mg1655/star_annotated',
    'BDGP6':
    '/home/cmb-panasas2/skchoudh/genomes/drosophila_melanogaster_BDGP6/star_annotated',
    'GRCz10':
    '/home/cmb-panasas2/skchoudh/genomes/GRCz10/star_annotated',
}

gtf_map = {
    'hg38':
    '/home/cmb-panasas2/skchoudh/genomes/hg38/annotation/gencode.v25.annotation.gtf',
    'mm10':
    '/home/cmb-panasas2/skchoudh/genomes/mm10/annotation/gencode.vM11.annotation.gtf',
    'sacCerR64':
    '/home/cmb-panasas2/skchoudh/genomes/sacCerR64/annotation/Saccharomyces_cerevisiae.R64-1-1.91.gtf',
    'MG1655':
    '/home/cmb-panasas2/skchoudh/genomes/escherichia_coli_str_k_12_substr_mg1655/annotation/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.38.gtf',
    'BDGP6':
    '/home/cmb-panasas2/skchoudh/genomes/drosophila_melanogaster_BDGP6/annotation/Drosophila_melanogaster.BDGP6.91.gtf',
    'GRCz10':
    '/home/cmb-panasas2/skchoudh/genomes/GRCz10/annotation/Danio_rerio.GRCz10.91.gtf'
}


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


def copy_sra_data(df,
                  sra_source_dir='/staging/as/skchoudh/SRA_datasets/',
                  sra_dest_dir='/staging/as/skchoudh/re-ribo-datasets/'):
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
        source_loc = srp_source_dir + os.path.sep + df_group[
            'experiment_accession'].str.cat(
                df_group['run_accession'] + '.sra', sep=os.path.sep)
        dest_loc = srp_dest_dir + os.path.sep + df_group[
            'experiment_accession'].str.cat(
                df_group['run_accession'] + '.sra', sep=os.path.sep)
        with tqdm(total=len(source_loc)) as pbar:
            for source, dest in zip(source_loc, dest_loc):
                mkdir_p(os.path.dirname(dest))
                if os.path.isfile(source):
                    symlink_force(source, dest)
                pbar.update()


def write_config(species, srp):
    rawdata_dir = os.path.join(samples_to_process_dir, species, srp)
    out_dir = os.path.join(re_ribo_analysis_dir, species, srp)
    gene_bed = os.path.join(riboraptor_annotation_dir, species,
                            genome_annotation_map[species], 'gene.bed.gz')
    utr5_bed = os.path.join(riboraptor_annotation_dir, species,
                            genome_annotation_map[species], 'utr5.bed.gz')
    utr3_bed = os.path.join(riboraptor_annotation_dir, species,
                            genome_annotation_map[species], 'utr3.bed.gz')
    intron_bed = os.path.join(riboraptor_annotation_dir, species,
                              genome_annotation_map[species], 'intron.bed.gz')
    start_codon_bed = os.path.join(riboraptor_annotation_dir, species,
                                   genome_annotation_map[species],
                                   'start_codon.bed.gz')
    stop_codon_bed = os.path.join(riboraptor_annotation_dir, species,
                                  genome_annotation_map[species],
                                  'stop_codon.bed.gz')
    cds_bed = os.path.join(riboraptor_annotation_dir, species,
                           genome_annotation_map[species], 'cds.bed.gz')
    genome_fasta = genome_fasta_map[species]
    gtf = gtf_map[species]
    chrom_sizes = chrom_sizes_map[species]
    star_index = star_index_map[species]
    to_write = """
    RAWDATA_DIR = '{}'
    OUT_DIR = '{}'
    GENOME_FASTA = '{}'
    CHROM_SIZES = '{}'
    STAR_INDEX = '{}'
    GTF = '{}'
    GENE_BED = '{}'
    STAR_CODON_BED = '{}'
    STOP_CODON_BED = '{}'
    CDS_BED = '{}'
    UTR5_BED = '{}'
    UTR3_BED = '{}'
    INTRON_BED = '{}'
    ORIENTATIONS = ['5prime', '3prime']
    STRANDS = ['pos', 'neg', 'combined']
    FRAGMENT_LENGTHS =  range(18, 39)
    """.format(rawdata_dir, out_dir, genome_fasta, chrom_sizes, star_index,
               gtf, gene_bed, start_codon_bed, stop_codon_bed, cds_bed,
               utr5_bed, utr3_bed, intron_bed)
    return dedent(to_write)


def create_config_file(df):
    df_grouped = df.groupby(['taxon_id'])

    for taxon_id, df_group in df_grouped:
        assert len(
            df_group['study_accession'].unique()) == 1, 'Multiple SRPs found'
        species = taxon_id_map[taxon_id]
        srp = df_group['study_accession'].unique()[0]
        with open(
                os.path.join(re_ribo_config_dir, '{}_{}.py'.format(
                    species, srp)), 'w') as fh:

            config = write_config(species, srp)
            fh.write(config)
            print('Wrote {}'.format(
                os.path.join(re_ribo_config_dir, '{}_{}.py'.format(
                    species, srp))))
