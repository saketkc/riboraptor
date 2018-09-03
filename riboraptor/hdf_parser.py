from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from collections import Counter
import os
import subprocess

import numpy as np
import pandas as pd
import pybedtools
import h5py
import six
import pyBigWig
from tqdm import tqdm
import warnings
import sys
from joblib import delayed

from .genome import _get_bed
from .genome import __GENOMES_DB__
from .interval import Interval

from .helpers import complementary_strand
from .helpers import mkdir_p
from .helpers import merge_intervals
from .helpers import path_leaf
from .parallel import ParallelExecutor
from .count import multiprocess_gene_coverage
import tempfile
TMP_DIR_ROOT = '/tmp'

#from multiprocessing import Pool


def _create_bigwig_from_bed(bed, chrom_lengths, outfile):
    """Create bigwig from a bed

    Parameters
    ----------
    bed: DataFrame
        Dataframe with standard six columns with appropriate names
    chrom_lengths: list of tuples
                   Sorted list of chromosome names and sizes
    outfile: string
             Path to store output file
    """
    if not isinstance(bed, pd.DataFrame):
        bed = pybedtools.BedTool(bed).as_dataframe()

    bw = pyBigWig.open(outfile, 'w')
    bw.addHeader(chrom_lengths, maxZooms=0)
    bw.addEntries(
        bed['chrom'].tolist(),
        bed['start'].tolist(),
        ends=bed['end'].tolist(),
        values=bed['score'].tolist())
    bw.close()


def hdf_to_bigwig(hdf, prefixdir):
    """Create fragment and strand specific bigwigs from hdf

    Parameters
    ----------
    hdf: string
         Path to hdf file
    prefix: string
            Prefix to store output files
    TODO: one fix for this would be to create separate nonoverlapping bigwigs
    and then merge them with bigWigCat (one more dependency)
    """
    hdf = h5py.File(hdf, 'r')
    chrom_names = list(map(lambda x: str(x), hdf['chrom_names']))
    protocol = hdf.attrs
    chrom_sizes = hdf['chrom_sizes']
    # This is already sorted in the file
    # So no need to sort again (as required by bigwig)
    chrom_lengths = list(zip(chrom_names, chrom_sizes))
    read_lengths = hdf['read_lengths']
    mkdir_p(prefixdir)
    for read_length in read_lengths:
        read_len_group = hdf['fragments'][read_length]
        for orientation in hdf['fragments'][read_length].keys():
            orientation_group = read_len_group[orientation]
            mkdir_p(os.path.join(prefixdir, str(read_length)))
            pos_bw = os.path.join(prefixdir, str(read_length), '{}_{}.bw'.format(orientation,
                                                                                 'pos'))
            neg_bw = os.path.join(prefixdir, str(read_length), '{}_{}.bw'.format(orientation,
                                                                                 'neg'))
            # This flle will store only the relevant
            # strand information
            collapsed_bw = os.path.join(prefixdir, str(read_length), '{}_{}.bw'.format(orientation,
                                                                                       'collapsed'))
            pos_bw = pyBigWig.open(pos_bw, 'w')
            neg_bw = pyBigWig.open(neg_bw, 'w')
            collapsed_bw = pyBigWig.open(collapsed_bw, 'w')
            pos_bw.addHeader(chrom_lengths, maxZooms=0)
            neg_bw.addHeader(chrom_lengths, maxZooms=0)
            collapsed_bw.addHeader(chrom_lengths, maxZooms=0)
            chrom_strand = list()
            for c in orientation_group.keys():
                chrom_strand.append((c.split('__')[0], c.split('__')[1]))
            chrom_strand = sorted(chrom_strand, key=lambda x: x[0])
            for chrom, mapped_strand in chrom_strand:

                chrom = str(chrom)
                chrom_strand = '{}__{}'.format(chrom, mapped_strand)
                assert mapped_strand in [
                    'pos', 'neg'
                ], 'Found strand {}'.format(mapped_strand)
                if mapped_strand == 'pos':
                    mapped_strand = '+'
                elif mapped_strand == 'neg':
                    mapped_strand = '-'

                chrom_obj = orientation_group[chrom_strand]
                starts = np.array(chrom_obj['positions'])
                ends = starts + 1
                values = np.array(chrom_obj['counts'])
                gene_strands = chrom_obj['gene_strand']
                if mapped_strand == '+':
                    pos_bw.addEntries(
                        [chrom] * len(starts),
                        starts,
                        ends=ends,
                        values=values)
                if mapped_strand == '+':
                    neg_bw.addEntries(
                        [chrom] * len(starts),
                        starts,
                        ends=ends,
                        values=values)
                gene_strand = chrom_obj['gene_strand']
                if protocol == 'forward':
                    # Should_keep mapped_strand + and gene_strand =
                    # OR mapped_strand - and gene_strand -i
                    for start, end, gene_strand, value in zip(
                            starts, ends, gene_strands, values):
                        if mapped_strand == gene_strand:
                            collapsed_bw.addEntries(
                                [chrom], [start], ends=[end], values=[value])

                elif protocol == 'reverse':
                    # Should keep reads mapped to + and gene_strand -
                    # OR mapped_strand - and gene_strand +
                    for start, end, gene_strand, value in zip(
                            starts, ends, gene_strands, values):
                        if mapped_strand == complementary_strand(gene_strand):
                            collapsed_bw.addEntries(
                                [chrom], [start], ends=[end], values=[value])

                elif protocol == 'unstranded':
                    # Count everything
                    # Do we still use this??
                    raise ValueError('Unstranded protocol not supported yet')

            pos_bw.close()
            neg_bw.close()
            collapsed_bw.close()
    hdf.close()


def tsv_to_bigwig(df, chrom_lengths, prefix):
    """Convert tsv (created by bam-coverage) to bigwig

    We create multiple bigwigs, each separately for fragment length,
    and the strand.

    This will output the following files, where N represents the fragment length

    N.5prime.pos.bw
    N.3prime.pos.bw

    N.5prime.neg.bw
    N.3prime.neg.bw


    Parameters
    ----------
    df: string or pd.DataFrame
        Path to bam-coverage outputted tsvi
    chrom_lengths: list of tuples
                   A list of tuples with chromsome name and size
                    [('chr1', 10000), ['chr2', 5000]]
    prefix: string
            Path to store the final tsv
    """
    if isinstance(chrom_lengths, six.string_types):
        if not os.path.isfile(chrom_lengths):
            raise RuntimeError(
                'Need either list of tuples or chrom.sizes file to proceed')
        sizes = []
        with open(chrom_lengths) as fh:
            for line in fh:
                chrom, size = line.strip().split('\t')
                sizes.append((str(chrom), int(size)))
        chrom_lengths = sizes
    chrom_lengths = list(sorted(chrom_lengths, key=lambda x: x[0]))
    if isinstance(df, six.string_types):
        df = pd.read_table(df)
    df = df.fillna(0)
    df.count_pos_strand = df.count_pos_strand.astype(float)
    df.count_neg_strand = df.count_neg_strand.astype(float)
    df_grouped = df.groupby('read_length')
    for read_length, read_len_group in df_grouped:
        orientation_grouped = read_len_group.groupby('orientation')
        for orientation, orientation_group in orientation_grouped:
            orientation_group = orientation_group.sort_values(
                by=['chrom', 'start', 'end']).drop(
                    columns=['read_length', 'orientation'])
            orientation_group_pos = orientation_group.copy().rename(
                columns={'count_pos_strand': 'score'})
            orientation_group_pos['strand'] = '+'
            orientation_group_pos = orientation_group_pos.sort_values(
                by=['chrom', 'start', 'end'])
            orientation_group_neg = orientation_group.copy().rename(
                columns={'count_neg_strand': 'score'})
            orientation_group_neg['strand'] = '-'
            orientation_group_neg = orientation_group_neg.sort_values(
                by=['chrom', 'start', 'end'])
            _create_bigwig_from_bed(
                orientation_group_pos, chrom_lengths, '{}_{}_{}_{}.bw'.format(
                    prefix, read_length, orientation, 'pos'))
            _create_bigwig_from_bed(
                orientation_group_neg, chrom_lengths, '{}_{}_{}_{}.bw'.format(
                    prefix, read_length, orientation, 'neg'))


def _process_gene_group(data):
    filepath, gene_group, chromosome_lengths, offset_5p, offset_3p, orientation, n_bases, fragment_length = data
    h5py_obj = HDFParser(filepath)
    regions = list(
        zip(gene_group['chrom'], gene_group['start'], gene_group['end'],
            gene_group['strand']))
    strand = regions[0][3]
    intervals = [Interval(i[0], i[1], i[2], i[3]) for i in regions]
    intervals_combined, gene_offset_5p, gene_offset_3p = merge_intervals(
        intervals, chromosome_lengths, offset_5p, offset_3p)
    regions = [(interval.chrom, interval.start, interval.end, interval.strand)
               for interval in intervals]
    all_dfs = [
        h5py_obj.get_coverage(region, fragment_length, orientation)
        for region in regions
    ]
    h5py_obj.close()
    all_dfs = pd.concat(all_dfs, ignore_index=True)
    # Reindex if the start position is '-'
    if strand == '-':
        all_dfs = all_dfs.reindex(np.arange(len(all_dfs.index) - 1, -1, -1))
    max_positions = min(len(all_dfs.index), n_bases)
    # Remove the 'start' column as we will
    # dill with index from now on
    coverage = all_dfs.loc[list(range(max_positions)), all_dfs.columns[1:]]
    return coverage


def _gene_group_to_tsv(data):
    gene_group = data['gene_group']
    offset_5p = data['offset_5p']
    offset_3p = data['offset_3p']
    #fragment_lengh = data['fragment_length']
    #orientation = data['orientation']
    chromosome_lengths = data['chromosome_lengths']
    n_bases = data['n_bases']
    tsv = data['tsv']
    tsv = tsv.set_index(['chrom', 'start', 'end'])
    regions = list(
        zip(gene_group['chrom'], gene_group['start'], gene_group['end'],
            gene_group['strand']))
    strand = regions[0][3]
    intervals = [Interval(i[0], i[1], i[2], i[3]) for i in regions]
    intervals_combined, gene_offset_5p, gene_offset_3p = merge_intervals(
        intervals, chromosome_lengths, offset_5p, offset_3p)
    regions = [(interval.chrom, interval.start, interval.end, interval.strand)
               for interval in intervals]
    df = []
    should_stop = False
    if strand == '+':
        total = 0
        for region in intervals_combined:
            if should_stop:
                break
            for position in range(region.start, region.end):
                if should_stop:
                    break
                #coverage_records =  tsv.loc[(region.chrom, position)]
                df.append(
                    [region.chrom, position, position + 1, region.strand])
                total += 1
                if total > n_bases:
                    should_stop = True
    else:
        # Iterate over the intervals in a reversed manner
        for region in intervals_combined[::-1]:
            if should_stop:
                break
            for position in np.arange(region.end - 1, region.end - 1, -1):
                if should_stop:
                    break
                #coverage_records =  tsv.loc[(region.chrom, position)]
                df.append(
                    [region.chrom, position, position + 1, region.strand])
                total += 1
                if total > n_bases:
                    should_stop = True

    df = pd.DataFrame(df)
    df.columns = ['chrom', 'start', 'end', 'strand']
    df = df.set_index(['chrom', 'start', 'end'])
    merged = df.join(tsv)
    merged = merged.reset_index()
    merged = merged.sort_values(by=['chrom', 'start', 'end'])
    if strand == '+':
        merged['index'] = range(len(merged.index))
    else:
        merged['index'] = np.range(len(merged.index), -1, -1)
    return merged


def create_metagene_from_multi_bigwig(bed,
                                      bigwigs,
                                      max_positions=1000,
                                      offset_5p=0,
                                      offset_3p=0,
                                      n_jobs=16,
                                      saveto=None):
    """Collapse multiple bigwigs to get bigwig.

    Test case -> we should be able to get the same coverage file as from export_metagene coverage
    if everything above this works okay.
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

    position_counter = Counter()
    metagene_coverage = pd.Series()
    for bw in tqdm(bigwigs):
        data = [(gene_group, bw.wig_location, offset_5p, offset_3p,
                 max_positions) for gene_name, gene_group in bed_grouped]
        aprun = ParallelExecutor(n_jobs=n_jobs)
        total = len(bed_grouped.groups)
        all_coverages = aprun(total=total)(
            delayed(multiprocess_gene_coverage)(d) for d in data)
        for norm_cov in all_coverages:
            metagene_coverage = metagene_coverage.add(norm_cov, fill_value=0)
            position_counter += Counter(norm_cov.index.tolist())

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


def merge_bigwigs(bigwigs, chrom_sizes, saveto, scale=False):
    """Merge multiple bigwigs into one.

    Note: This seems impossible doing it through pybigWig way
    and hence the dependency on ucsc-bigWigMerge

    Parameters
    ----------
    bigwigs: list
             List of path to bigwigs
    chrom_size: string
                Path to chromsome.sizes file
    saveto: string
            Path for outputting bigwig

    """
    assert len(bigwigs) > 1, 'Need more bigwigs to merge'
    with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
        bedgraph = os.path.join(temp_dir, '{}.bedGraph'.format(
            path_leaf(saveto)))
        cmds = ['bigWigMerge'] + bigwigs + [bedgraph]
        try:
            p = subprocess.Popen(
                cmds,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True)
            stdout, stderr = p.communicate()
        except FileNotFoundError:
            raise FileNotFoundError(
                "bigWigMerge not found on the path. This is an external "
                "tool from UCSC which can be downloaded from "
                "http://hgdownload.soe.ucsc.edu/admin/exe/. Alternatatively, use "
                "`conda install ucsc-bigwigmerge`")
        cmds = ['bedSort', bedgraph, '{}.sorted'.format(bedgraph)]
        try:
            p = subprocess.Popen(
                cmds,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True)
            stdout, stderr = p.communicate()
        except FileNotFoundError:
            raise FileNotFoundError(
                "bedSort not found on the path. This is an external "
                "tool from UCSC which can be downloaded from "
                "http://hgdownload.soe.ucsc.edu/admin/exe/. Alternatatively, use "
                "`conda install ucsc-bedsort`")

        cmds = [
            'bedGraphToBigWig', '{}.sorted'.format(bedgraph), chrom_sizes,
            saveto
        ]
        try:
            p = subprocess.Popen(
                cmds,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True)
            stdout, stderr = p.communicate()
        except FileNotFoundError:
            raise FileNotFoundError(
                "bedGraphToBigwig not found on the path. This is an external "
                "tool from UCSC which can be downloaded from "
                "http://hgdownload.soe.ucsc.edu/admin/exe/. Alternatatively, use "
                "`conda install ucsc-bedgraphtobigwig`")


class HDFParser(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.h5py_obj = h5py.File(filepath, 'r')
        chrom_names = self.h5py_obj['chrom_names']
        chrom_sizes = self.h5py_obj['chrom_sizes']
        chrom_lengths = dict(
            [(chrom.replace('_neg', '').replace('_pos', ''), int(size))
             for chrom, size in zip(chrom_names, chrom_sizes)])
        self.chromosome_lengths = chrom_lengths

    def close(self):
        self.h5py_obj.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.h5py_obj.close()

    def get_coverage(self,
                     region=(None, None, None, None),
                     fragment_length=None,
                     orientation='5prime',
                     outprefix=None):
        """Get coverage for a selected region.

        Since the region is user provided, we assume the user is aware if it is strand
        specfic or not

        Paramters
        ---------
        region: tuple
                (chromosome, start[0-based], end[1-based], strand(+/-))
        fragment_length: list or string
                         list of fragment lengths to use
                         Example: [20, 21, 22] or 'all'

        """
        chrom, start, stop, strand = region
        start = int(start)
        stop = int(stop)
        assert chrom is not None, 'chromosome not set'
        assert start is not None, 'start position is not set'
        assert stop is not None, 'end position not set'
        assert strand is not None, 'strand is not set'
        assert start < stop, 'start should be < stop'
        assert fragment_length is not None, 'fragment_length not set'
        if isinstance(fragment_length, int):
            fragment_length = [fragment_length]

        if fragment_length == 'all':
            h5py_fragments = self.h5py_obj['read_lengths']
            fragment_length = h5py_fragments
        fragment_length = list(map(lambda x: str(x), fragment_length))
        coverages_normalized = pd.DataFrame()
        coverages = pd.DataFrame()
        position_counter = Counter()
        for l in fragment_length:
            if strand == '-':
                chrom_name = '{}__neg'.format(chrom)
            elif strand == '+':
                chrom_name = '{}__pos'.format(chrom)
            else:
                raise ValueError('strand ill-defined: {}'.format(strand))
            root_obj = self.h5py_obj['fragments'][l][orientation]
            if chrom_name not in root_obj.keys():
                # This chromosome does not exist in the
                # key value store
                # So should returns zeros all way
                coverage = pd.Series(
                    [0] * (stop - start), index=range(start, stop))

            else:
                chrom_obj = root_obj[chrom_name]
                counts_series = pd.Series(
                    list(chrom_obj['counts']),
                    index=list(chrom_obj['positions']))
                #print(counts_series)
                coverage = counts_series.get(list(range(start, stop)))
                if coverage is None:
                    coverage = pd.Series(
                        [0] * (stop - start), index=range(start, stop))
                coverage = coverage.fillna(0)
            # Mean is taken by summing the rows
            coverage_mean = coverage.mean(axis=0, skipna=True)
            # to normalize
            # we divide the sum by vector obtained in previous
            # such that each column gets divided
            coverage_normalized = coverage.divide(coverage_mean).fillna(0)
            coverages = coverages.join(
                pd.DataFrame(coverage, columns=[str(l)]), how='outer')
            coverages_normalized = coverages_normalized.join(
                pd.DataFrame(coverage_normalized, columns=[str(l)]),
                how='outer')
            position_counter += Counter(coverage.index.tolist())
        position_counter = pd.Series(Counter(position_counter)).sort_index()
        coverage_sum = coverages.sum(axis=1)
        coverage_normalized_sum = coverages_normalized.sum(axis=1)
        coverage_sum_div = coverage_sum.divide(position_counter, axis='index')
        coverage_normalized_sum_div = coverage_normalized_sum.divide(
            position_counter, axis='index')
        if outprefix:
            mkdir_p(os.path.dirname(outprefix))
            coverage_sum_div.to_csv(
                '{}_raw.tsv', index=False, header=True, sep='\t')
            coverage_normalized_sum_div.to_csv(
                '{}_normalized.tsv', index=False, header=True, sep='\t')

        coverages.index.name = 'start'
        coverages = coverages.reset_index()
        coverages_normalized.index.name = 'start'
        coverages_normalized = coverages_normalized.reset_index()
        return coverages, coverages_normalized, coverage_sum_div, coverage_normalized_sum_div

    def get_read_length_dist(self):
        read_lengths = list(
            map(lambda x: int(x), self.h5py_obj['read_lengths']))
        read_counts = list(self.h5py_obj['read_lengths_counts'])
        return pd.Series(read_counts, index=read_lengths).sort_index()
