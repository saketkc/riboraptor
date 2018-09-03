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
import subprocess
import sys

import numpy as np
import pandas as pd
import pybedtools
import pysam
import six
import h5py
from tqdm import tqdm

from .genome import _get_sizes
from .genome import _get_bed
from .genome import __GENOMES_DB__

from .helpers import merge_intervals
from .helpers import mkdir_p
from .helpers import complementary_strand

from .wig import WigReader
from .interval import Interval
from .parallel import ParallelExecutor
from joblib import delayed
from .infer_protocol import infer_protocol
from .helpers import read_bed_as_intervaltree
from .helpers import is_read_uniq_mapping, create_bam_index
from .helpers import find_first_non_none, find_last_non_none


class OrderedCounter(Counter, OrderedDict):
    'Counter that remembers the order elements are first encountered'

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, OrderedDict(self))

    def __reduce__(self):
        return self.__class__, (OrderedDict(self), )


def _get_gene_strand(refseq, chrom, mapped_start, mapped_end):
    """Lookup a particular location on chromosome to determine
    its strand


    Parameters
    ----------
    refseq: string or intervaltree
            refseq genes.bed file loaded as

    """
    if isinstance(refseq, six.string_types):
        refseq = read_bed_as_intervaltree(refseq)
    gene_strand = list(set(refseq[chrom].find(mapped_start, mapped_end)))
    if len(gene_strand) > 1:
        return 'ambiguous'
    if len(gene_strand) == 0:
        # Try searching upstream 30 and downstream 30 nt
        gene_strand = list(
            set(refseq[chrom].find(mapped_start - 30, mapped_end - 30)))
        if len(gene_strand) == 0:
            # Downstream
            gene_strand = list(
                set(refseq[chrom].find(mapped_start + 30, mapped_end + 30)))
            if len(gene_strand) == 0:
                return 'not_found'
            if len(gene_strand) > 1:
                return 'ambiguous'
            return gene_strand[0]
        if len(gene_strand) > 1:
            return 'ambiguous'
        return gene_strand[0]
    return gene_strand[0]


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
    original_interval_length: int
                              Total length of the entire gene interval (not accounting for any offsets)
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

    # if it is located on negative strand
    # reverse the values since we no longer
    # want to track the individual position
    # but just an indexed version
    if strand == '-':
        coverages_combined.reverse()
    coverages_combined = np.array(coverages_combined).flatten()
    original_interval_length = len(coverages_combined) - gene_offset_5p - gene_offset_3p
    coverages_combined = pd.Series(
        coverages_combined,
        index=np.arange(-gene_offset_5p,
                        len(coverages_combined) - gene_offset_5p))

    return (coverages_combined, gene_offset_5p, gene_offset_3p, original_interval_length)


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


def multiprocess_gene_coverage(data):
    """Process gene_c overage given a bigwig and a genegroup.

    WigReader is not pickleable when passed as an argument so we use strings
    as input for the bigwig

    Parameters
    ----------
    data: tuple
          gene_gorup, bigwig, offset_5p, offset_3p, max_positions, orientation

    Returns
    -------
    norm_cov: Series
              normalized coverage
    """
    gene_group, bw, offset_5p, offset_3p, max_positions, orientation = data
    bw = WigReader(bw)
    coverage, gene_offset_5p, gene_offset_3p, original_gene_length = gene_coverage(gene_group, bw, offset_5p, offset_3p)
    coverage = coverage.fillna(0)


    if orientation == '5prime':
        if max_positions is not None and len(coverage.index) > 0:
            # min_index will correspond to the gene_offset_5p in general
            min_index = min(coverage.index.tolist())
            max_index = max(coverage.index.tolist())
            assert min_index == -gene_offset_5p, 'min_index and gene_offset_5p are not same| min_index: {} | gene_offset_5p: {}'.format(min_index,
                                                                                                                                        -gene_offset_5p)
            coverage = coverage.get(np.arange(min_index, min(max_index,
                                                        max_positions)))
    elif orientation == '3prime':
        # We now want to be tracking things from the end position
        # we can do this since gene_coverage() takes care of the strand
        # so a 3prime is always the tail of the array
        # note that if gene_offset_5p >0, in this case, it is almost never used
        # since we restrict ourselves to max_positions, which itself is almost
        # always < 1000
        if max_positions is not None and len(coverage.index) > 0:
            max_index = max(coverage.index.tolist())
            min_index = min(coverage.index.tolist())
            assert min_index == -gene_offset_5p, 'min_index and gene_offset_5p are not same| min_index: {} | gene_offset_5p: {}'.format(min_index,
                                                                                                                                        -gene_offset_5p)
            # max_index is the maximum we can go to the right
            # our stop codon will be located gene_offset_3p upstream of this index
            # Let's reindex our series so that we set
            coverage = coverage.reindex(np.arange(-max_index, -min_index, 1))
            coverage = coverage.get(np.arange(-max_positions, gene_offset_3p))
    else:
        raise ValueError('{} orientation not supported'.format(orientation))

    assert coverage is not None, 'coverage is none | max_index={} | min_index={}| gene_offset_3p={} | gene_offset_5p={}'.format(max_index, min_index,
                                                                                                                                gene_offset_3p, gene_offset_5p)
    coverage_mean = coverage.mean()
    norm_cov = coverage / coverage_mean
    norm_cov = norm_cov.fillna(0)
    bw.close()
    return norm_cov


def export_metagene_coverage(bed,
                             bw,
                             max_positions=None,
                             saveto=None,
                             offset_5p=0,
                             offset_3p=0,
                             orientation='5prime',
                             n_jobs=16):
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
    orientation: string
                 ['5prime', '3prime'] indicating the end of read
                 being tracked
    n_jobs: int
            Number of paralle threads to open
            Better to do on a multi-cpu machine, but also works decently on
            a single core machine

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
    data = [(gene_group, bw.wig_location, offset_5p, offset_3p, max_positions, orientation)
            for gene_name, gene_group in bed_grouped]
    aprun = ParallelExecutor(n_jobs=n_jobs)
    total = len(bed_grouped.groups)
    all_coverages = aprun(total=total)(
        delayed(multiprocess_gene_coverage)(d) for d in data)
    for norm_cov in all_coverages:
        metagene_coverage = metagene_coverage.add(norm_cov, fill_value=0)
        position_counter += Counter(norm_cov.index.tolist())
    del all_coverages
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
    create_bam_index(bam)
    bam = pysam.AlignmentFile(bam, 'rb')
    read_counts = Counter([
        read.query_length for read in bam.fetch() if is_read_uniq_mapping(read)
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
    create_bam_index(bam)
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
    create_bam_index(bam)
    bam = pysam.AlignmentFile(bam, 'rb')
    n_mapped = 0
    for read in bam.fetch():
        if is_read_uniq_mapping(read):
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
    create_bam_index(inbam)
    allreadsbam = pysam.AlignmentFile(inbam, 'rb')
    uniquereadsbam = pysam.AlignmentFile(outbam, 'wb', template=allreadsbam)
    total_count = allreadsbam.count()
    with tqdm(total=total_count) as pbar:
        for read in allreadsbam.fetch():
            if is_read_uniq_mapping(read):
                uniquereadsbam.write(read)
            pbar.update()
    allreadsbam.close()
    uniquereadsbam.close()


def get_bam_coverage(bam, bed, outprefix=None):
    """Get coverage from bam

    Parameters
    ----------
    bam: string
         Path to bam file
    bed: string
         Path to genes.bed or cds.bed file required for inferring protocol of bam

    Returns
    -------
    coverage: dict(dict)
              dict with keys as chrom:position with query lengths as keys
    """
    protocol, forward_mapped_reads, reverse_mapped_reads, total_reads = infer_protocol(
        bam, bed)
    if isinstance(bam, six.string_types):
        create_bam_index(bam)
        bam = pysam.AlignmentFile(bam, 'rb')

    coverage = defaultdict(lambda: defaultdict(OrderedCounter))
    total_counts = bam.count()
    mismatches = defaultdict(lambda: defaultdict(OrderedCounter))
    with tqdm(total=total_counts) as pbar:
        for read in bam.fetch():
            if not is_read_uniq_mapping(read):
                pbar.update()
                continue
            if read.is_reverse:
                strand = '-'
            else:
                strand = '+'
            # Why full_length?
            # If full_length is set, None values will be included for any soft-clipped or unaligned positions within the read.i
            # The returned list will thus be of the same length as the read.
            reference_pos = read.get_reference_positions(full_length=True)
            # Store 5' position
            position_5prime = None
            # Store 3' position
            position_3prime = None

            # We track 0-based start
            # the start coordinates are always 0-based
            # irrepective of the strand

            ######################################################
            """
            Forward strand:

            -------------------------
            0                       24
                  -------------------
                  5'                3'
            5' position = 6 (0-based)
            3' position = 24 (0-based)
            Read length = 24-6+1 = 19

            Reverse strand
            -------------------------
            0                       24
                  -------------------
                  3'                5'

            5' position = 24 (0-based)
            3' position = 6 (0-based)

            """
            # query_alignment_start = read.query_alignment_start
            # query_alignment_end = read.query_alignment_end

            query_alignment_length = read.query_alignment_length

            mismatches_or_softclipping = False
            query_length = read.query_length
            if query_alignment_length < query_length:
                # the aligned length of this
                # read is shorter than the original
                # read length
                # could happen because of soft-clipping or mismatches
                # STAR and most other RNA-seq mappers do not recommend using a masked
                # genome and we don't too, so these will almost always
                # be mismatches
                mismatches_or_softclipping = True
            assert query_length == len(
                reference_pos
            ), 'reference_pos and query_length should be equal'
            mismatch_at_3prime = False
            mismatch_at_5prime = False
            if strand == '+':
                position_5prime = reference_pos[0]
                position_3prime = reference_pos[-1]
                if position_5prime is None:
                    # print('query_alignment_start: {} | read_length: {}'.format(query_alignment_start, query_length))
                    # Mismatch/sofclip at 5'
                    # Need to do some adjusting
                    # query_alignment_start refers to the first position on the read whihch has a non-mismatch non-softclipped
                    # base
                    # This is a bit ahead in the read so we subtract
                    # position_5prime = reference_pos[query_alignment_end] - query_alignment_start
                    idx, position = find_first_non_none(reference_pos)
                    position_5prime = position - idx
                    mismatch_at_5prime = True
                    if position_5prime <0:
                        sys.stderr.write('position_5prime<0 | idx : {} | position: {} | reference_pos: {}'.format(idx, position, reference_pos))
                        sys.exit(1)
                if position_3prime is None:
                    # print('query_alignment_end: {} | read_length: {}'.format(query_alignment_end, query_length))
                    # position_3prime = reference_pos[query_alignment_end] + query_length - query_alignment_end - 1

                    # We cannot do the folowing since this would force no splicing
                    # position_3prime = position_5prime + query_length - 1
                    idx, position = find_last_non_none(reference_pos)
                    position_3prime = position + idx
                    mismatch_at_3prime = True
                #   assert position_3prime == position_5prime + query_length - 1, 'Position prime not equal? {} vs {}'.format(
                #   position_3prime, position_5prime + query_length - 1)

            else:
                # Negative strand so no need to adjust
                # switch things
                position_5prime = reference_pos[-1]
                position_3prime = reference_pos[0]

                if position_5prime is None:
                    # print(reference_pos, read.reference_end)
                    # print('query_alignment_end: {} | read_length: {}'.format(query_alignment_end, query_length))
                    # position_5prime = reference_pos[query_alignment_end] + query_length - query_alignment_end - 1
                    idx, position = find_last_non_none(reference_pos)
                    position_5prime = position + idx
                    mismatch_at_5prime = True
                    # print(position_5prime)
                    #sys.exit(1)
                if position_3prime is None:
                    # print('query_alignment_start: {} | read_length: {}'.format(query_alignment_start, query_length))
                    # position_3prime = reference_pos[query_alignment_start] - query_alignment_start
                    idx, position = find_first_non_none(reference_pos)
                    position_3prime = position - idx
                    mismatch_at_3prime = True
                    if position_3prime <0:
                        sys.stderr.write('position_3prime<0 | idx : {} | position: {} | reference_pos: {}'.format(idx, position, reference_pos))
                        sys.exit(1)
                    # print(position_3prime)
                    #sys.exit(1)
                #   assert position_5prime == position_3prime + query_length - 1, 'Position prime not equal? {} vs {}'.format(
                #    position_5prime, position_3prime + query_length - 1)

            assert position_5prime >= 0, 'Wrong 5prime position: {}'.format(
                position_5prime)
            assert position_3prime >= 0, 'Wrong 3prime position: {}'.format(
                position_3prime)

            coverage['{}:{}:5prime'.format(
                read.reference_name, position_5prime)][query_length]['+'] += 0
            coverage['{}:{}:5prime'.format(
                read.reference_name, position_5prime)][query_length]['-'] += 0
            coverage['{}:{}:3prime'.format(
                read.reference_name, position_3prime)][query_length]['+'] += 0
            coverage['{}:{}:3prime'.format(
                read.reference_name, position_3prime)][query_length]['-'] += 0

            coverage['{}:{}:5prime'.format(
                read.reference_name,
                position_5prime)][query_length][strand] += 1
            coverage['{}:{}:3prime'.format(
                read.reference_name,
                position_3prime)][query_length][strand] += 1

            mismatches['{}:{}:3prime:{}'.format(
                read.reference_name, position_3prime,
                strand)][query_length]['match'] += int(not mismatch_at_3prime)
            mismatches['{}:{}:3prime:{}'.format(
                read.reference_name, position_3prime,
                strand)][query_length]['mismatch'] += int(mismatch_at_3prime)

            mismatches['{}:{}:5prime:{}'.format(
                read.reference_name, position_5prime,
                strand)][query_length]['match'] += int(not mismatch_at_5prime)
            mismatches['{}:{}:5prime:{}'.format(
                read.reference_name, position_5prime,
                strand)][query_length]['mismatch'] += int(mismatch_at_5prime)

            pbar.update()
    if outprefix:
        df = pd.DataFrame.from_dict(
            {(i, j): coverage[i][j]
             for i in coverage.keys() for j in coverage[i].keys()},
            orient='index')
        df = df.reset_index()

        mismatches_df = pd.DataFrame.from_dict(
            {(i, j): mismatches[i][j]
             for i in mismatches.keys() for j in mismatches[i].keys()},
            orient='index')
        mismatches_df = mismatches_df.reset_index()
        mismatches_df.columns = [
            'chr_pos_orient_strand', 'read_length', 'match', 'mismatch'
        ]

        mismatches_df[['chrom', 'start', 'orientation', 'strand'
                       ]] = mismatches_df['chr_pos_orient_strand'].str.split(
                           ':', n=-1, expand=True)
        #assert list(df.copy().columns)[::-1][0:2] == ['-', '+'], 'column orders out of order: {}'.format(list(df.columns)[::-1])
        mismatches_df = mismatches_df.drop(columns=['chr_pos_orient_strand'])
        mismatches_df = mismatches_df[[
            'chrom', 'start', 'orientation', 'strand', 'read_length', 'match',
            'mismatch'
        ]]
        mismatches_df['start'] = mismatches_df['start'].astype(int)
        mismatches_df = mismatches_df.sort_values(
            by=['chrom', 'start', 'read_length', 'orientation'])
        mismatches_df.to_csv(
            '{}_mismatches.tsv'.format(outprefix),
            sep='\t',
            index=False,
            header=True)
        """
        Stored as:

            chrom\tstart_position(0-based)\tnumber of hits on + strand\tnumber of hits on - strand
        """
        df.columns = [
            'chr_pos_orient',
            'read_length',
            'count_pos_strand',
            'count_neg_strand',
        ]
        # n=-1 tto expand all splits
        df[['chrom', 'start', 'orientation']] = df['chr_pos_orient'].str.split(
            ':', n=-1, expand=True)
        df['start'] = df['start'].astype(int)
        df['count_pos_strand'] = df['count_pos_strand'].fillna(0).astype(int)
        df['count_neg_strand'] = df['count_neg_strand'].fillna(0).astype(int)
        df['end'] = df['start'] + 1
        df = df[[
            'chrom', 'start', 'end', 'read_length', 'orientation',
            'count_pos_strand', 'count_neg_strand'
        ]]
        bed_tree = read_bed_as_intervaltree(bed)
        df['gene_strand'] = df.apply(
            lambda row: _get_gene_strand(bed_tree, row['chrom'], row['start'], row['end']),
            axis=1)
        df = df.sort_values(by=[
            'chrom', 'start', 'end', 'read_length', 'orientation',
            'count_pos_strand', 'count_neg_strand', 'gene_strand'
        ])

        df.to_csv(
            '{}.tsv'.format(outprefix), sep='\t', index=False, header=True)

        reference_and_length = OrderedDict(
            sorted(
                zip(bam.header.references, bam.header.lengths),
                key=lambda x: x[0]))

        df = df.sort_values(by=[
            'read_length', 'chrom', 'start', 'orientation', 'count_pos_strand'
        ])
        df = df.rename(columns={
            'count_pos_strand': 'pos',
            'count_neg_strand': 'neg'
        })
        df_molten = pd.melt(
            df,
            id_vars=[
                'chrom', 'start', 'end', 'read_length', 'orientation',
                'gene_strand'
            ],
            value_vars=['pos', 'neg'])
        df_molten = df_molten.rename(columns={'variable': 'strand'})
        df_molten = df_molten.sort_values(by=[
            'chrom', 'start', 'end', 'read_length', 'orientation', 'strand',
            'gene_strand'
        ])
        df_molten['chrom'] = df_molten.chrom.str.cat(
            df_molten.strand, sep='__')
        # Convert read length dtype to make it more compatible to h5py (requires string keys)
        df_molten['read_length'] = df_molten.read_length.astype(str)
        df_molten = df_molten.drop(columns=['end', 'strand'])
        df_readlen_grouped = df_molten.groupby(['read_length'])

        df_readlen_orient = pd.DataFrame(
            df_molten.groupby(['read_length',
                               'orientation']).value.agg('sum')).reset_index()
        # Check if the counts are same irrespective of orientation
        df_readlen_orient_unmelt = df_readlen_orient.pivot(
            index='read_length', columns='orientation', values='value')
        assert (df_readlen_orient_unmelt['5prime'] == df_readlen_orient_unmelt[
            '3prime']).all()
        dt = h5py.special_dtype(vlen=str)
        read_lengths_list = df_readlen_orient_unmelt.index.tolist()
        read_lengths_counts = df_readlen_orient_unmelt.loc[read_lengths_list,
                                                           '5prime'].tolist()

        with tqdm(total=len(df_readlen_grouped)) as pbar:
            with h5py.File('{}.hdf5'.format(outprefix), 'w') as h5py_file:
                h5py_file.attrs['protocol'] = protocol
                dset = h5py_file.create_dataset(
                    'read_lengths', (len(read_lengths_list), ),
                    dtype=dt,
                    compression='gzip',
                    compression_opts=9)
                dset[...] = np.array(read_lengths_list)

                dset = h5py_file.create_dataset(
                    'read_lengths_counts', (len(read_lengths_list), ),
                    dtype=np.dtype('int64'),
                    compression='gzip',
                    compression_opts=9)
                dset[...] = np.array(read_lengths_counts)

                dset = h5py_file.create_dataset(
                    'chrom_names', (len(reference_and_length.keys()), ),
                    dtype=dt,
                    compression='gzip',
                    compression_opts=9)
                dset[...] = np.array(list(reference_and_length.keys()))

                dset = h5py_file.create_dataset(
                    'chrom_sizes', (len(reference_and_length.values()), ),
                    dtype=np.dtype('int64'),
                    compression='gzip',
                    compression_opts=9)
                dset[...] = np.array(list(reference_and_length.values()))

                h5_readlen_group = h5py_file.create_group('fragments')

                for read_length, df_readlen_group in df_readlen_grouped:
                    h5_readlen_subgroup = h5_readlen_group.create_group(
                        read_length)
                    df_orientation_grouped = df_readlen_group.groupby(
                        'orientation')
                    for orientation, df_orientation_group in df_orientation_grouped:
                        h5_orientation_subgroup = h5_readlen_subgroup.create_group(
                            orientation)
                        df_chrom_grouped = df_orientation_group.groupby(
                            'chrom')
                        for chrom, chrom_group in df_chrom_grouped:
                            counts_series = pd.Series(
                                chrom_group['value'].tolist())
                            positions_series = pd.Series(
                                chrom_group['start'].tolist())
                            chrom_group = h5_orientation_subgroup.create_group(
                                chrom)
                            dset = chrom_group.create_dataset(
                                'positions', (len(counts_series), ),
                                dtype=np.dtype('int64'),
                                compression='gzip',
                                compression_opts=9)
                            dset[...] = positions_series

                            dset = chrom_group.create_dataset(
                                'counts', (len(counts_series), ),
                                dtype=np.dtype('float64'),
                                compression='gzip',
                                compression_opts=9)
                            dset[...] = counts_series

                            dset = chrom_group.create_dataset(
                                'gene_strand', (len(counts_series), ),
                                dtype=dt,
                                compression='gzip',
                                compression_opts=9)
                            dset[...] = chrom_group['gene_strand']
                    pbar.update()
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
