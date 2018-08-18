"""Utilities for read counting operations.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings

from collections import defaultdict
from collections import Counter
from collections import OrderedDict
from functools import reduce
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd
import pybedtools
import pysam
import six
from tqdm import tqdm

from .genome import _get_sizes
from .genome import _get_bed
from .genome import __GENOMES_DB__

from .helpers import merge_intervals
from .helpers import mkdir_p
from .helpers import complementary_strand

from .wig import WigReader
from .interval import Interval

# Unmapped, Unmapped+Reverse strand, Not primary alignment,
# Not primary alignment + reverse strand, supplementary alignment

# Source: https://broadinstitute.github.io/picard/explain-flags.html
__SAM_NOT_UNIQ_FLAGS__ = [4, 20, 256, 272, 2048]


def _create_bam_index(bam):
    """Create bam index.

    Parameters
    ----------
    bam : str
          Path to bam file
    """
    if not os.path.exists('{}.bai'.format(bam)):
        pysam.index(bam)


def _is_read_uniq_mapping(read):
    """Check if read is uniquely mappable.

    Parameters
    ----------
    read : pysam.Alignment.fetch object


    Most reliable: ['NH'] tag
    """
    # Filter out secondary alignments
    if read.is_secondary:
        return False
    tags = dict(read.get_tags())
    try:
        nh_count = tags['NH']
    except KeyError:
        # Reliable in case of STAR
        if read.mapping_quality == 255:
            return True
        if read.mapping_quality < 1:
            return False
        # NH tag not set so rely on flags
        if read.flag in __SAM_NOT_UNIQ_FLAGS__:
            return False
        else:
            raise RuntimeError('Malformed BAM?')
    if nh_count == 1:
        return True
    return False


def gene_coverage(gene_group, bw, offset_5p=0, offset_3p=0):
    """Get gene coverage.

    Parameters
    ----------
    gene_group: DataFrame
                gene group from bed file
    bw: str
        Path to bigwig to fetch the scores from
    offset_5p: int (positive)
               Number of bases to count upstream (5')
    offset_3p: int (positive)
               Number of bases to count downstream (3')

    Returns
    -------
    coverage_combined: series
                       Series with index as position and value as coverage
    gene_offset_5p: Gene wise 5 prime offset
                    This might be different from `offset_5p` in cases where
                    `offset_5p` leads to a negative coordinate
    gene_offset_3p: Gene wise 3 prime offset
                    This might be different from `offset_3p` in cases where
                    `offset_3p` leads to position beyond chromsome length
    """
    if offset_5p < 0 or offset_3p < 0:
        raise RuntimeError('Offsets must be non-negative')
        sys.exit(1)
    if not isinstance(bw, WigReader):
        bw = WigReader(bw)
    chromosome_lengths = bw.chromosomes

    if len(gene_group['strand'].unique()) != 1:
        raise RuntimeError('Multiple strands?: {}'.format(gene_group))
    if len(gene_group['chrom'].unique()) != 1:
        raise RuntimeError('Chrom not unique for: {}'.format(gene_group))
    strand = gene_group['strand'].unique()[0]

    intervals = list(
        zip(gene_group['chrom'], gene_group['start'], gene_group['end'],
            gene_group['strand']))

    intervals = [Interval(i[0], i[1], i[2], i[3]) for i in intervals]
    intervals_combined, gene_offset_5p, gene_offset_3p = merge_intervals(
        intervals, chromosome_lengths, offset_5p, offset_3p)
    coverages = bw.query(intervals_combined)
    if len(coverages) == 0:
        return (pd.Series([]), 0, 0)
    coverages_combined = []
    for cov in coverages:
        coverages_combined += list(cov)
    if strand == '-':
        coverages_combined.reverse()
    coverages_combined = np.array(coverages_combined).flatten()
    coverages_combined = pd.Series(
        coverages_combined,
        index=np.arange(-gene_offset_5p,
                        len(coverages_combined) - gene_offset_5p))
    return (coverages_combined, gene_offset_5p, gene_offset_3p)


def export_gene_coverages(bed, bw, saveto, offset_5p=0, offset_3p=0):
    """Export all gene coverages.

    Parameters
    ----------
    bed: str
         Path to CDS or 5'UTR or 3'UTR bed
    bw: str
        Path to bigwig to fetch the scores from
    saveto: str
            Path to write output tsv file
    offset_5p: int (positive)
               Number of bases to count upstream (5')
    offset_3p: int (positive)
               Number of bases to count downstream (3')

    Returns
    -------
    gene_profiles: file
                   with the following format:
                   gene1\t5poffset1\t3poffset1\tcnt1_1 cnt1_2 cnt1_3 ...\n
                   gene2\t5poffset2\t3poffset2\tcnt2_1 cnt2_2 cnt2_3 ...\n
    """
    if bed.lower().split('_')[0] in __GENOMES_DB__:
        splitted = bed.lower().split('_')
        if len(splitted) == 2:
            genome, region_type = splitted
        elif len(splitted) == 3:
            genome = splitted[0]
            region_type = ('_').join(splitted[1:])
        bed = _get_bed(region_type, genome)
    bed_df = pybedtools.BedTool(bed).sort().to_dataframe()
    bed_df['chrom'] = bed_df['chrom'].astype(str)
    bed_df['name'] = bed_df['name'].astype(str)
    bed_grouped = bed_df.groupby('name')

    if not isinstance(bw, WigReader):
        bw = WigReader(bw)

    to_write = 'gene_name\toffset_5p\toffset_3p\tcoverage\n'
    for gene_name, gene_group in tqdm(bed_grouped):
        coverage, gene_offset_5p, gene_offset_3p = gene_coverage(
            gene_group, bw, offset_5p, offset_3p)
        coverage = coverage.fillna(0)
        coverage = coverage.astype(int)
        coverage = coverage.tolist()

        to_write += '{}\t{}\t{}\t{}\n'.format(gene_name, int(gene_offset_5p),
                                              int(gene_offset_3p), coverage)

    mkdir_p(os.path.dirname(saveto))
    with open(saveto, 'w') as outfile:
        outfile.write(to_write)


def export_metagene_coverage(bed,
                             bw,
                             max_positions=None,
                             saveto=None,
                             offset_5p=0,
                             offset_3p=0):
    """Export metagene coverage.

    Parameters
    ----------
    bed: str
         Path to CDS or 5'UTR or 3'UTR bed
    bw: str
        Path to bigwig to fetch the scores from
    max_positions: int
                   Number of positions to consider while
                   calculating the normalized coverage
                   Higher values lead to slower implementation
    saveto: str
            Path to write output tsv file
    offset_5p: int (positive)
               Number of bases to count upstream (5')
    offset_3p: int (positive)
               Number of bases to count downstream (3')

    Returns
    -------
    metagene_coverage: series
                       Metagene coverage
    """
    if max_positions is None:
        warnings.warn("The max_positions is not specified, it could"
                      "take long time to calculate metagene coverage")
    if max_positions is not None and max_positions <= 0:
        raise RuntimeError('The max_positions must be positive')
        sys.exit(1)

    if bed.lower().split('_')[0] in __GENOMES_DB__:
        splitted = bed.lower().split('_')
        if len(splitted) == 2:
            genome, region_type = splitted
        elif len(splitted) == 3:
            genome = splitted[0]
            region_type = ('_').join(splitted[1:])
        bed = _get_bed(region_type, genome)
    bed_df = pybedtools.BedTool(bed).sort().to_dataframe()
    bed_df['chrom'] = bed_df['chrom'].astype(str)
    bed_df['name'] = bed_df['name'].astype(str)
    bed_grouped = bed_df.groupby('name')

    if not isinstance(bw, WigReader):
        bw = WigReader(bw)

    position_counter = Counter()
    metagene_coverage = pd.Series()

    for gene_name, gene_group in tqdm(bed_grouped):
        coverage, _, _ = gene_coverage(gene_group, bw, offset_5p, offset_3p)
        coverage = coverage.fillna(0)

        if max_positions is not None and len(coverage.index) > 0:
            min_index = min(coverage.index.tolist())
            max_index = max(coverage.index.tolist())
            coverage = coverage[np.arange(min_index,
                                          min(max_index, max_positions))]
        coverage_mean = coverage.mean()
        if coverage_mean > 0:
            norm_cov = coverage / coverage_mean
            metagene_coverage = metagene_coverage.add(norm_cov, fill_value=0)
            position_counter += Counter(coverage.index.tolist())

    if len(position_counter) != len(metagene_coverage):
        raise RuntimeError('Gene normalizaed counter mismatch')
        sys.exit(1)

    position_counter = pd.Series(position_counter)
    metagene_coverage = metagene_coverage.div(position_counter)

    if saveto:
        mkdir_p(os.path.dirname(saveto))
        to_write = pd.DataFrame({
            'position': metagene_coverage.index,
            'count': metagene_coverage.values
        })
        to_write = to_write[['position', 'count']]
        to_write.to_csv(saveto, sep=str('\t'), index=False)

    return metagene_coverage


def export_read_counts(gene_coverages, saveto, keep_offsets=True):
    """export read counts from gene coverages file.

    Parameters
    ----------
    gene_coverages: string
                    Path to gene coverages.tsv
    saveto: str
            Path to save output tsv
            gene_name\tcount\tlength
    keep_offsets: bool
                 whether to keep the 5' and 3' offsets in gene coverages
                 default is False
    """
    if '.gz' in gene_coverages:
        gene_coverages_df = pd.read_table(gene_coverages, compression='gzip')
    else:
        gene_coverages_df = pd.read_table(gene_coverages)
    gene_coverages_zipped = list(
        zip(gene_coverages_df['gene_name'], gene_coverages_df['offset_5p'],
            gene_coverages_df['offset_3p'], gene_coverages_df['coverage']))
    to_write = "gene_name\tcount\tlength\n"
    for gene_name, offset_5p, offset_3p, cov in gene_coverages_zipped:
        coverage = eval(cov)
        coverage = pd.Series(
            np.array(coverage),
            index=np.arange(-offset_5p,
                            len(coverage) - offset_5p))
        coverage = coverage.fillna(0)
        max_index = max(coverage.index.tolist())
        if not keep_offsets:
            coverage = coverage[np.arange(0, max_index - offset_3p)]
        count = coverage.sum()
        length = len(coverage.tolist())
        to_write += "{}\t{}\t{}\n".format(gene_name, int(count), int(length))

    mkdir_p(os.path.dirname(saveto))
    with open(saveto, 'w') as output:
        output.write(to_write)


def merge_gene_coverages(gene_coverages, max_positions=None, saveto=None):
    """merge gene coverages to generate metagene coverage.

    Parameters
    ----------
    gene_coverages: string
                    Path to gene coverages.tsv
    max_positions: int
                   Number of positions to consider while
                   calculating the normalized coverage
    saveto: str
            Path to save output tsv

    Returns
    -------
    metagene_coverage : Series
                        Metagene coverage
    """
    if '.gz' in gene_coverages:
        gene_coverages_df = pd.read_table(gene_coverages, compression='gzip')
    else:
        gene_coverages_df = pd.read_table(gene_coverages)
    gene_coverages_zipped = list(
        zip(gene_coverages_df['offset_5p'], gene_coverages_df['offset_3p'],
            gene_coverages_df['coverage']))

    position_counter = Counter()
    metagene_coverage = pd.Series()
    for offset_5p, offset_3p, cov in gene_coverages_zipped:
        coverage = eval(cov)
        coverage = pd.Series(
            np.array(coverage),
            index=np.arange(-offset_5p,
                            len(coverage) - offset_5p))
        coverage = coverage.fillna(0)

        if max_positions is not None and len(coverage.index) > 0:
            min_index = min(coverage.index.tolist())
            max_index = max(coverage.index.tolist())
            coverage = coverage[np.arange(min_index,
                                          min(max_index, max_positions))]
        coverage_mean = coverage.mean()
        if coverage_mean > 0:
            norm_cov = coverage / coverage_mean
            metagene_coverage = metagene_coverage.add(norm_cov, fill_value=0)
            position_counter += Counter(coverage.index.tolist())

    if len(position_counter) != len(metagene_coverage):
        raise RuntimeError('Gene normalized counter mismatch')
        sys.exit(1)

    position_counter = pd.Series(position_counter)
    metagene_coverage = metagene_coverage.div(position_counter)

    if saveto:
        mkdir_p(os.path.dirname(saveto))
        to_write = pd.DataFrame({
            'position': metagene_coverage.index,
            'count': metagene_coverage.values
        })
        to_write = to_write[['position', 'count']]
        to_write.to_csv(saveto, sep=str('\t'), index=False)

    return metagene_coverage


def merge_read_counts(read_counts, saveto):
    """merge read counts tsv files to one count table

    Parameters
    ----------
    read_counts: str
                 Path to file which contains paths of all read counts
                 tsv file you want to merge
                 format:
                 sample1 path1
                 sample2 path2
    saveto: str
            Path to save output tsv

    """
    with open(read_counts) as f:
        samples = f.readlines()
    samples = [x.strip().split() for x in samples]
    dfs = []
    for sample, path in samples:
        df = pd.read_table(path, index_col=0)
        df = df[['count']]
        df.rename(columns={'count': sample}, inplace=True)
        dfs.append(df)
    df_merged = reduce(lambda left, right: pd.merge(
                                             left, right, how='outer',
                                             left_index=True,
                                             right_index=True), dfs)
    df_merged.fillna(0, inplace=True)
    df_merged = df_merged.astype(int)
    df_merged.to_csv(saveto, sep=str('\t'))


def export_read_length(bam, saveto=None):
    """Count read lengths.

    Parameters
    ----------
    bam: str
         Path to bam file
    saveto: str
            Path to write output tsv file
    Returns
    -------
    lengths: counter
             Counter of read length and counts

    """
    _create_bam_index(bam)
    bam = pysam.AlignmentFile(bam, 'rb')
    read_counts = Counter([
        read.query_length for read in bam.fetch()
        if _is_read_uniq_mapping(read)
    ])
    if saveto:
        mkdir_p(os.path.dirname(saveto))
        to_write = 'read_length\tcount\n'
        for read_length, count in six.iteritems(read_counts):
            to_write += '{}\t{}\n'.format(read_length, count)
        with open(saveto, 'w') as output:
            output.write(to_write)
    return read_counts


def read_enrichment(read_lengths, min_length=28, max_length=32):
    """Calculate read enrichment for a certain range of lengths

    Parameters
    ----------
    read_lengths: str
                  Path to read length tsv
    min_length: int
                The low end of the range
    max_length: int
                The high end of the range

    Returns
    -------
    ratio: float
           Enrichment in this range (Scale 0-1)

    """

    if '.gz' in read_lengths:
        read_lengths_df = pd.read_table(
            read_lengths,
            names=['frag_len', 'frag_cnt'],
            sep='\t',
            compression='gzip')
    else:
        read_lengths_df = pd.read_table(
            read_lengths, names=['frag_len', 'frag_cnt'], sep='\t')
    read_lengths = pd.Series(
        read_lengths_df.frag_cnt.tolist(),
        index=read_lengths_df.frag_len.tolist())

    target = read_lengths[np.arange(min_length, max_length + 1)].sum()
    total = read_lengths.sum()
    ratio = target / total
    return ratio


def get_region_sizes(bed):
    """Get collapsed lengths of gene in bed.

    Parameters
    ----------
    bed: str
         Path to bed file

    Returns
    -------
    region_sizes: dict
                  Region sizes with gene names as key
                  and value as size of this named region
    """
    if bed.lower().split('_')[0] in __GENOMES_DB__:
        genome, region_type = bed.lower().split('_')
        bed = _get_bed(region_type, genome)
    bed = pybedtools.BedTool(bed).to_dataframe()
    bed['chrom'] = bed['chrom'].astype(str)
    bed['name'] = bed['name'].astype(str)
    bed_grouped = bed.groupby('name')
    region_sizes = {}
    for gene_name, gene_group in bed_grouped:
        intervals = list(
            zip(gene_group['chrom'], gene_group['start'], gene_group['end'],
                gene_group['strand']))
        intervals = [Interval(i[0], i[1], i[2], i[3]) for i in intervals]
        intervals_combined, _, _ = merge_intervals(intervals)
        for interval in intervals_combined:
            if gene_name not in region_sizes:
                region_sizes[gene_name] = interval.end - interval.start
            else:
                region_sizes[gene_name] += interval.end - interval.start

    return region_sizes


def bedgraph_to_bigwig(bedgraph, sizes, saveto, input_is_stream=False):
    """Convert bedgraph to bigwig.

    Parameters
    ----------
    bedgraph : str
               Path to bedgraph file
    sizes : str
            Path to genome chromosome sizes file
            or genome name
    saveto : str
             Path to write bigwig file
    input_is_stream : bool
                      True if input is sent through stdin
    """
    if input_is_stream:
        total_lines = len(bedgraph)
        with open(os.path.splitext(saveto)[0] + '.bg', 'w') as fp:
            for index, line in enumerate(bedgraph):
                if index == (total_lines - 1):
                    fp.write(line.rstrip())
                else:
                    fp.write(line)
            filename = str(fp.name)
        bedgraph = filename

    if not os.path.isfile(sizes):
        if sizes in __GENOMES_DB__:
            sizes = _get_sizes(sizes)
        else:
            raise RuntimeError('Could not load size for {}'.format(sizes))

    cmds = ['bedSort', bedgraph, bedgraph]
    p = subprocess.Popen(
        cmds,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True)
    stdout, stderr = p.communicate()
    rc = p.returncode
    if rc != 0:
        raise RuntimeError(
            'Error running bedSort.\nstdout : {} \n stderr : {}'.format(
                stdout, stderr))

    cmds = ['bedGraphToBigWig', bedgraph, sizes, saveto]
    p = subprocess.Popen(
        cmds,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True)
    stdout, stderr = p.communicate()
    rc = p.returncode
    if rc != 0:
        raise RuntimeError(
            'Error running bedSort.\nstdout : {} \n stderr : {}'.format(
                stdout, stderr))


def bam_to_bedgraph(bam, strand='both', end_type='5prime', saveto=None):
    """Create bigwig from bam.

    Parameters
    ----------
    bam : str
          Path to bam file
    strand : str, optional
             Use reads mapping to '+/-/both' strands
    end_type : str
               Use only end_type=5prime(5') or "3prime(3')"
    saveto : str, optional
              Path to write bedgraph

    Returns
    -------
    genome_cov : str
                 Bedgraph output

    """
    if strand not in ['+', '-', 'both']:
        raise RuntimeError('Strand should be one of \'+\', \'-\', \'both\'')
    if end_type == '5prime':
        extra_args = '-5'
    elif end_type == '3prime':
        extra_args = '-3'
    elif end_type == 'either':
        extra_args = ''
    bed = pybedtools.BedTool(bam)
    if strand != 'both':
        genome_cov = bed.genome_coverage(
            bg=True, strand=strand, additional_args=extra_args)
    else:
        genome_cov = bed.genome_coverage(bg=True, additional_args=extra_args)
    if saveto:
        with open(saveto, 'w') as outf:
            outf.write(str(genome_cov))
    return str(genome_cov)


def mapping_reads_summary(bam, saveto=None):
    """Count number of mapped reads.

    Parameters
    ----------
    bam: str
         Path to bam file
    saveto: str
            Path to write output tsv

    Returns
    -------
    counts : counter
             Counter with keys as number of times read maps
             and values as number of reads of that type
    """
    _create_bam_index(bam)
    bam = pysam.AlignmentFile(bam, 'rb')
    counts = Counter()
    for read in bam.fetch():
        if read.is_secondary:
            continue
        try:
            nh_count = Counter([dict(read.get_tags())['NH']])
        except KeyError:
            nh_count = Counter([1])
        counts += nh_count

    if saveto:
        mkdir_p(os.path.dirname(saveto))
        to_write = 'n_mapping_times\tn_reads\n'
        for n_mapped, n_reads in six.iteritems(counts):
            to_write += '{}\t{}\n'.format(n_mapped, n_reads)
        with open(saveto, 'w') as output:
            output.write(to_write)
    return counts


def count_uniq_mapping_reads(bam):
    """Count number of mapped reads.

    Parameters
    ----------
    bam: str
         Path to bam file

    Returns
    -------
    n_mapped: int
              Count of uniquely mapped reads
    """
    _create_bam_index(bam)
    bam = pysam.AlignmentFile(bam, 'rb')
    n_mapped = 0
    for read in bam.fetch():
        if _is_read_uniq_mapping(read):
            n_mapped += 1
    bam.close()
    return n_mapped


def extract_uniq_mapping_reads(inbam, outbam):
    """Extract only uniquely mapping reads from a bam.

    Parameters
    ----------
    inbam: string
           Path to input bam file
    outbam: string
            Path to write unique reads bam to
    """
    _create_bam_index(inbam)
    allreadsbam = pysam.AlignmentFile(inbam, 'rb')
    uniquereadsbam = pysam.AlignmentFile(outbam, 'wb', template=allreadsbam)
    total_count = allreadsbam.count()
    with tqdm(total=total_count) as pbar:
        for read in allreadsbam.fetch():
            if _is_read_uniq_mapping(read):
                uniquereadsbam.write(read)
            pbar.update()
    allreadsbam.close()
    uniquereadsbam.close()


def get_bam_coverage(bam,
                     orientation='5prime',
                     saveto=None):
    """ Get coverage from bam given orientation

    Parameters
    ----------
    bam: string
         Path to bam file

    orientation: string
                 5prime/3prime/both

    Returns
    -------
    coverage: dict(dict)
              dict with keys as chrom:position with query lengths as keys
    """
    if isinstance(bam, six.string_types):
        bam = pysam.AlignmentFile(bam, 'rb')

    coverage = defaultdict(lambda: defaultdict(Counter))
    total_counts = bam.count()
    with tqdm(total=total_counts) as pbar:
        for read in bam.fetch():
            if not _is_read_uniq_mapping(read):
                pbar.update()
                continue
            if read.is_reverse:
                strand = '-'
            else:
                strand = '+'
            reference_pos = read.get_reference_positions()
            if strand == '+':
                if orientation == '5prime':
                    # Track 5' end
                    position = reference_pos[0]
                elif orientation == '3prime':
                    # Track 3' end
                    position = reference_pos[-1]

            else:
                # Negative strand so no need to adjust
                # switch things
                if orientation == '5prime':
                    # Track 5' end on negative strand
                    position = reference_pos[-1]
                elif orientation == '3prime':
                    # Track 3' end on negative strand
                    position = reference_pos[0]
            query_length = read.query_length
            coverage['{}:{}'.format(read.reference_name,
                                    position)][query_length][strand] += 1
            pbar.update()
    if saveto:
        df = pd.DataFrame.from_dict(
            {(i, j): coverage[i][j]
             for i in coverage.keys() for j in coverage[i].keys()},
            orient='index')
        df = df.reset_index()
        """
        Stored as:
            chrom\tstart_position(0-based)\tnumber of hits on + strand\tnumber of hits on - strand
        """
        df.columns = [
            'chr_pos', 'read_length', 'count_pos_strand', 'count_neg_strand'
        ]
        df[['chrom', 'start']] = df['chr_pos'].str.split(':', n=1, expand=True)
        df['start'] = df['start'].astype(int)
        df['count_pos_strand'] = df['count_pos_strand'].fillna(0).astype(int)
        df['count_neg_strand'] = df['count_neg_strand'].fillna(0).astype(int)
        df['end'] = df['start'] + 1
        df = df[[
            'chrom', 'start', 'end', 'read_length', 'count_pos_strand',
            'count_neg_strand'
        ]]
        df = df.sort_values(
            by=['chrom', 'start', 'read_length', 'count_pos_strand'])
        df.to_csv(saveto, sep='\t', index=False, header=True)
        return df
    return coverage


def get_bam_coverage_on_bed(bam,
                            bed,
                            protocol='forward',
                            orientation='5prime',
                            max_positions=1000,
                            offset=60,
                            saveto=None):
    """Get bam coverage over start_codon/stop_codon coordinates

    Parameter
    ---------
    bam: string
         Path to bam file
    bed: string
         Path to bed file
    protocol: string
              forward/unstranded/reverse
              can be gussed from infer-protocol
    orientation: string
                 5prime/3prime ends to be tracked
    max_positions: int
                   Number of positions to consider for metagene counting
                   Larger is slower(!)
    offset: int
            Positive if to be count upstream of 5prime
            Does not support 3prime stuff properly yet
    saveto: string
            path to store the output tsv
    """

    assert protocol in ['forward', 'reverse', 'unstranded']
    if bed.lower().split('_')[0] in __GENOMES_DB__:
        splitted = bed.lower().split('_')
        if len(splitted) == 2:
            genome, region_type = splitted
        elif len(splitted) == 3:
            genome = splitted[0]
            region_type = ('_').join(splitted[1:])
        bed = _get_bed(region_type, genome)
    bed_df = pybedtools.BedTool(bed).sort().to_dataframe()
    bed_df['chrom'] = bed_df['chrom'].astype(str)
    bed_df['name'] = bed_df['name'].astype(str)
    pos_names = bed_df.chrom.str.cat(bed_df.start.astype(str), sep=':')
    pos_names = pos_names.str.cat(bed_df.strand.astype(str), sep=':').tolist()

    print('Fetching coverage:')
    coverage_dict = get_bam_coverage(bam, orientation=orientation, saveto=None)
    metagene_coverage = defaultdict(Counter)
    print('Collapsing to metagene:')
    for pos_name in tqdm(pos_names):
        chrom, start_position, gene_strand = pos_name.split(':')
        start_position = int(start_position)
        cur_pointer = -offset
        while cur_pointer < max_positions:
            curr_pos = start_position + cur_pointer
            counts_counter = coverage_dict['{}:{}'.format(chrom, curr_pos)]
            query_counter = OrderedDict()

            # Only keep relevant strand
            for query_length, strand_dict in six.iteritems(counts_counter):
                if protocol == 'unstranded':
                    # sum both strand mapping reads
                    query_counter[query_length] = sum(strand_dict.values())
                elif protocol == 'forward':
                    # only keep reads mapping to same strand as gene strand
                    query_counter[query_length] = strand_dict[gene_strand]
                elif protocol == 'reverse':
                    # only keep reads mapping to opposite strand
                    query_counter[query_length] = strand_dict[
                        complementary_strand(gene_strand)]
                else:
                    raise ValueError(
                        'Not a valid protocol: {}'.format(protocol))
            metagene_coverage[cur_pointer] += query_counter
            cur_pointer += 1
    df = pd.DataFrame.from_dict(metagene_coverage)
    df = df.fillna(0)
    df = df.astype(int)
    if saveto:
        df.to_csv(saveto, sep='\t', index=True, header=True)
    return df

