from __future__ import absolute_import
from __future__ import print_function
import pytest
import numpy as np

from riboraptor.interval import Interval
from riboraptor.helpers import merge_intervals
from riboraptor.helpers import read_refseq_bed
from riboraptor.helpers import scale_bigwig
from riboraptor.hdf_parser import normalize_bw_hdf
from riboraptor.wig import WigReader


def test_scale_bigwig():
    bw = 'tests/data/SRX2536403_subsampled.unique.bigWig'
    scaled_bw = 'tests/data/SRX2536403_subsampled.unique.scaled.bigWig'
    chrom_sizes = 'riboraptor/annotation/hg38/hg38.chrom.sizes'
    scale_bigwig(bw, chrom_sizes, scaled_bw, 0.1)
    bw1 = WigReader(bw)
    bw2 = WigReader(scaled_bw)
    region = [Interval('chrY', 10197344, 10197345)]
    np.testing.assert_array_almost_equal(
        bw1.query(region) * 0.1, bw2.query(region))


def test_merge_intervals():
    bw = WigReader('tests/data/SRX2536403_subsampled.unique.bigWig')
    chromosome_lengths = bw.chromosomes

    # test when intervals are empty
    intervals = []
    intervals_combined = merge_intervals(intervals, chromosome_lengths, 30, 30)
    assert (intervals_combined == ([], 30, 30))

    # test when only one interval
    intervals = [Interval('chr1', 10, 20, '-')]
    intervals_combined = merge_intervals(intervals, chromosome_lengths, 30, 30)
    assert (intervals_combined == ([Interval('chr1', 0, 50, '-')], 30, 10))

    # test overlapping intervals
    intervals = [
        Interval('chr1', x[0], x[1], '-')
        for x in [(1, 9), (7, 15), (30, 40), (13, 18), (3, 4)]
    ]
    intervals_combined = merge_intervals(intervals, chromosome_lengths, 30, 30)
    expected = [Interval('chr1', x[0], x[1], '-') for x in [(0, 18), (30, 70)]]
    assert (intervals_combined == (expected, 30, 1))

    # test more tricky intervals
    intervals = [
        Interval('chr1', x[0], x[1], '-')
        for x in [(260, 400), (200, 300), (280, 300)]
    ]
    intervals_combined = merge_intervals(intervals, chromosome_lengths, 30, 30)
    assert (intervals_combined == ([Interval('chr1', 170, 430, '-')], 30, 30))

    # test intervals near the end of chrom for neg strand
    chr1_length = chromosome_lengths['chr1']
    intervals = [
        Interval('chr1', x[0], x[1], '-')
        for x in [(chr1_length - 10,
                   chr1_length - 2), (260, 400), (200, 300), (280, 300)]
    ]
    intervals_combined = merge_intervals(intervals, chromosome_lengths, 30, 30)
    assert (intervals_combined == ([
        Interval('chr1', 170, 400, '-'),
        Interval('chr1', chr1_length - 10, chr1_length, '-')
    ], 2, 30))

    # test intervals near the end of chrom for pos strand
    chr1_length = chromosome_lengths['chr1']
    intervals = [
        Interval('chr1', x[0], x[1], '+')
        for x in [(chr1_length - 10,
                   chr1_length - 2), (260, 400), (200, 300), (280, 300)]
    ]
    intervals_combined = merge_intervals(intervals, chromosome_lengths, 10, 30)
    assert (intervals_combined == ([
        Interval('chr1', 190, 400, '+'),
        Interval('chr1', chr1_length - 10, chr1_length, '+')
    ], 10, 2))

    # test intervals of one-based near start
    chr1_length = chromosome_lengths['chr1']
    intervals = [
        Interval('chr1', x[0], x[1], '+')
        for x in [(chr1_length - 10,
                   chr1_length - 2), (260, 400), (5, 300), (280, 300)]
    ]
    intervals_combined = merge_intervals(intervals, chromosome_lengths, 10, 30,
                                         False)
    assert (intervals_combined == ([
        Interval('chr1', 1, 400, '+'),
        Interval('chr1', chr1_length - 10, chr1_length, '+')
    ], 4, 2))


def test_refseq_read():
    refseq = read_refseq_bed('tests/data/hg38_v24_refseq.bed12')
    assert list(sorted(refseq.keys())) == list(
        sorted(['chr1', 'chr22_KI270734v1_random', 'chr22_KI270733v1_random']))
