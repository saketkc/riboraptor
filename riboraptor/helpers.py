"""All functions that are not so useful, but still useful."""

from collections import Counter
from collections import OrderedDict
from collections import defaultdict
import errno
import itertools
import math
import os
import re
import sys
import ntpath
import pickle
import subprocess

from scipy import stats
import numpy as np
import pandas as pd
import six
import pybedtools
import pysam
import pyBigWig
from bx.intervals.intersection import IntervalTree

import warnings
from .interval import Interval

# Unmapped, Unmapped+Reverse strand, Not primary alignment,
# Not primary alignment + reverse strand, supplementary alignment

# Source: https://broadinstitute.github.io/picard/explain-flags.html
__SAM_NOT_UNIQ_FLAGS__ = [4, 20, 256, 272, 2048]

CBB_PALETTE = [
    "#000000",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
]


def order_dataframe(df, columns):
    """Order a dataframe

    Order a dataframe by moving the `columns` in the front

    Parameters
    ----------
    df: Dataframe
        Dataframe
    columns: list
             List of columns that need to be put in front
    """
    if isinstance(columns, six.string_types):
        columns = [columns]  # let the command take a string or list
    remaining_columns = [w for w in df.columns if w not in columns]
    df = df[columns + remaining_columns]
    return df


def _fix_bed_coltype(bed):
    """Fix bed chrom and name columns to be string

    This is necessary since the chromosome numbers are often interpreted as int
    """
    bed["chrom"] = bed["chrom"].astype(str)
    bed["name"] = bed["name"].astype(str)
    return bed


def check_file_exists(filepath):
    """Check if file exists.

    Parameters
    ----------
    filepath : str
               Path to file
    """
    if os.path.isfile(os.path.abspath(filepath)):
        return True
    return False


def list_to_ranges(list_of_int):
    """Convert a list to a list of range object

    Parameters
    ----------

    list_of_int: list
        List of integers to be squeezed into range

    Returns
    -------

    list_of_range: list
        List of range objects


    """
    sorted_list = sorted(set(list_of_int))
    for key, group in itertools.groupby(enumerate(sorted_list), lambda x: x[1] - x[0]):
        group = list(group)
        yield group[0][1], group[-1][1]


def create_ideal_periodic_signal(signal_length):
    """Create ideal ribo-seq signal.

    Parameters
    ----------
    signal_length : int
                    Length of signal to create

    Returns
    -------
    signal : array_like
             1-0-0 signal

    """
    uniform_signal = np.array([4 / 6.0] * signal_length)
    uniform_signal[list(range(1, len(uniform_signal), 3))] = 1 / 6.0
    uniform_signal[list(range(2, len(uniform_signal), 3))] = 1 / 6.0
    return uniform_signal


def identify_peaks(coverage):
    """Given coverage array, find the site of maximum density"""
    return np.argmax(coverage[list(range(-18, -10))])


def millify(n):
    """Convert integer to human readable format.

    Parameters
    ----------
    n : int

    Returns
    -------
    millidx : str
              Formatted integer
    """
    if n is None or np.isnan(n):
        return "NaN"
    millnames = ["", " K", " M", " B", " T"]
    # Source: http://stackoverflow.com/a/3155023/756986
    n = float(n)
    millidx = max(
        0,
        min(
            len(millnames) - 1, int(math.floor(0 if n == 0 else math.log10(abs(n)) / 3))
        ),
    )

    return "{:.1f}{}".format(n / 10 ** (3 * millidx), millnames[millidx])


def mkdir_p(path):
    """Python version mkdir -p

    Parameters
    ----------

    path : str
    """
    if path:
        try:
            os.makedirs(path)
        except OSError as exc:  # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise


def symlink_force(source, destination):
    """Create forcelink forcefully

    Parameters
    ----------
    source: string
            Location to source file
    destination: string
                 Location to target

    """
    try:
        os.symlink(source, destination)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            os.remove(destination)
            os.symlink(source, destination)
        else:
            raise exc


def r2(x, y):
    """Calculate pearson correlation between two vectors.

    Parameters
    ----------
    x : array_like
        Input
    y : array_like
        Input
    """
    return stats.pearsonr(x, y)[0] ** 2


def round_to_nearest(x, base=5):
    """Round to nearest base.

    Parameters
    ----------
    x : float
        Input

    Returns
    -------
    v : int
        Output
    """
    return int(base * round(float(x) / base))


def set_xrotation(ax, degrees):
    """Rotate labels on x-axis.

    Parameters
    ----------
    ax : matplotlib.Axes
         Axes object
    degrees : int
              Rotation degrees
    """
    for i in ax.get_xticklabels():
        i.set_rotation(degrees)


def summary_stats_two_arrays_welch(
    old_mean_array,
    new_array,
    old_var_array=None,
    old_n_counter=None,
    carried_forward_observations=None,
):
    """Average two arrays using welch's method

    Parameters
    ----------
    old_mean_array : Series
                Series of previous means with index as positions
    old_var_array : Series
                Series of previous variances with index as positions
    new_array : array like
                Series of new observations
                (Does noes
                Ciunts of number of positions at a certain index

    Returns
    -------
    m : array like
        Column wise Mean array
    var : array like
         Column wise variance

    Consider an example: [1,2,3], [1,2,3,4], [1,2,3,4,5]

    old = [1,2,3]
    new = [1,2,3,4]
    counter = [1,1,1]
    mean = [1,2,3,4] Var =[na, na, na, na], carried_fowrad = [[1,1], [2,2], [3,3], [4]]

    old = [1,2,3,4]
    new = [1,2,3,4,5]
    couter = [2,2,2,1]

    mean = [1,2,3,4,5]
    var = [0,0,0, na, na]
    carried_forward = [[], [], [], [4,4], [5]]
    """
    if not isinstance(old_mean_array, pd.Series):
        old_mean_array = pd.Series(old_mean_array)
    if not isinstance(new_array, pd.Series):
        new_array = pd.Series(new_array)
    if old_n_counter is not None and not isinstance(old_n_counter, pd.Series):
        old_n_counter = pd.Series(old_n_counter)

    len_old, len_new = len(old_mean_array), len(new_array)
    if old_n_counter is None:
        # Initlaized from current series
        old_n_counter = pd.Series(
            np.zeros(len(old_mean_array)) + 1, index=old_mean_array.index
        )
    if old_var_array is None:
        # Initlaized from current series
        old_var_array = pd.Series(
            np.zeros(len(old_mean_array)) + np.nan, index=old_mean_array.index
        )
    # Update positions counts based on new_array
    new_n_counter = old_n_counter.add(
        pd.Series(np.zeros(len(new_array)) + 1, index=new_array.index), fill_value=0
    )
    if len_old > len_new:
        len_diff = len_old - len_new
        # Pad the incoming array
        # We append NAs to the end of new_array since it will mostly be in the metagene context
        max_index = np.max(new_array.index.tolist())
        new_index = np.arange(max_index + 1, max_index + 1 + len_diff)
        new_array = new_array.append(
            pd.Series(np.zeros(len_diff) + np.nan, index=new_index),
            verify_integrity=True,
        )
    elif len_old < len_new:
        len_diff = len_new - len_old
        # Pad the old array
        if len_old == 0:
            old_mean_array = pd.Series([])
        else:
            max_index = np.max(old_mean_array.index.tolist())
            new_index = np.arange(max_index + 1, max_index + 1 + len_diff)
            old_mean_array = old_mean_array.append(
                pd.Series(np.zeros(len_diff) + np.nan, index=new_index),
                verify_integrity=True,
            )

    if not (old_mean_array.index == new_array.index).all():
        print("old array index: {}".format(old_mean_array))
        print("new array index: {}".format(new_array))
    positions_with_less_than3_obs = defaultdict(list)
    for index, counts in six.iteritems(new_n_counter):
        # Which positions has <3 counts for calculating variance
        if counts <= 3:
            # Fetch the exact observations from history
            try:
                last_observations = carried_forward_observations[index]
            except:
                # No carreid forward passed
                if not np.isnan(old_mean_array[index]):
                    last_observations = [old_mean_array[index]]
                else:
                    last_observations = []
            # Add entry from new_array only if it is not NAN
            if not np.isnan(new_array[index]):
                last_observations.append(new_array[index])
            positions_with_less_than3_obs[index] = last_observations

    # positions_with_less_than3_obs = pd.Series(positions_with_less_than3_obs)

    # delta = x_n - mean(x_{n-1})
    delta = new_array.subtract(old_mean_array)
    """
    for index, value in six.iteritems( delta ):
        if np.isnan(value):
            if not np.isnan(old_mean_array[index]):
                delta[index] = old_mean_array[index]
            else:
                delta[index] = new_array[index]
    """

    # delta = delta/n
    delta_normalized = delta.divide(new_n_counter)
    # mean(x_n) = mean(x_{n-1}) + delta/n
    new_mean_array = old_mean_array.add(delta_normalized)
    for index, value in six.iteritems(new_mean_array):
        if np.isnan(value):
            if not np.isnan(old_mean_array[index]):
                new_mean_array[index] = old_mean_array[index]
            else:
                new_mean_array[index] = new_array[index]
    # print(delta)
    # print(new_n_counter)
    # print(delta_normalized)
    # print(new_mean_array)
    # mean_difference_current = x_n - mean(x_n)
    # mean_difference_previous = x_n - mean(x_{n-1})

    mean_difference_current = new_array.fillna(0) - new_mean_array.fillna(0)
    mean_difference_previous = new_array.fillna(0) - old_mean_array.fillna(0)

    # (x_n-mean(x_n))(x_n-mean(x_{n-1})
    product = np.multiply(mean_difference_current, mean_difference_previous)

    # (n-1)S_n^2 - (n-2)S_{n-1}^2 = (x_n-mean(x_n)) (x_n-mean(x_{n-1}))

    # old_ssq = (n-1)S_{n-1}^2
    # (n-2)S_{n-1}^2
    old_sum_of_sq = (old_n_counter - 2).multiply(old_var_array.fillna(0))

    # new_ssq = (old_ssq + product)
    # (n-1) S_n^2
    new_sum_of_sq = old_sum_of_sq + product

    # if counts is less than 3, set sum of sq to NA
    new_sum_of_sq[new_n_counter < 3] = np.nan

    # if counts just became 3, compute the variance
    for index, counts in six.iteritems(new_n_counter):
        if counts == 3:
            observations = positions_with_less_than3_obs[index]
            variance = np.var(observations)
            print(index, variance)
            new_sum_of_sq[index] = variance
            # delete it from the history
            del positions_with_less_than3_obs[index]

    new_var_array = new_sum_of_sq.divide(new_n_counter - 1)
    new_var_array[new_var_array == np.inf] = np.nan
    new_var_array[new_n_counter < 3] = np.nan
    """
    for index, counts in six.iteritems(new_n_counter):
        if counts < 3:
            if not np.isnan(new_array[index]):
                if index not in list(positions_with_less_than3_obs.keys()):
                    positions_with_less_than3_obs[index] = list()
                assert index in positions_with_less_than3_obs.keys()
                positions_with_less_than3_obs[index].append(new_array[index])
    """
    return new_mean_array, new_var_array, new_n_counter, positions_with_less_than3_obs


def path_leaf(path):
    """Get path's tail from a filepath"""
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def parse_star_logs(infile, outfile=None):
    """Parse star logs into a dict

    Parameters
    ----------
    infile : str
        Path to starlogs.final.out file

    Returns
    -------
    star_info : dict
                Dict with necessary records parsed
    """
    ANNOTATIONS = [
        "total_reads",
        "uniquely_mapped",
        "uniquely_mapped_percent",
        "multi_mapped_percent",
        "unmapped_percent",
        "multi_mapped",
    ]
    star_info = OrderedDict()
    with open(infile) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("Number of input reads"):
                star_info[ANNOTATIONS[0]] = int(line.strip().split("\t")[1])
            elif line.startswith("Uniquely mapped reads number"):
                star_info[ANNOTATIONS[1]] = int(line.strip().split("\t")[1])
            elif line.startswith("Uniquely mapped reads %"):
                star_info[ANNOTATIONS[2]] = round(
                    float(line.strip("%").split("\t")[1]), 2
                )
            elif line.startswith("Number of reads mapped to multiple loci"):
                star_info[ANNOTATIONS[5]] = int(line.strip().split("\t")[1])
            elif line.startswith("Number of reads mapped to too many loci"):
                star_info[ANNOTATIONS[5]] += int(line.strip().split("\t")[1])
            elif line.startswith("% of reads mapped to multiple loci"):
                star_info[ANNOTATIONS[3]] = round(
                    float(line.strip("%").split("\t")[1]), 2
                )
            elif line.startswith("% of reads mapped to too many loci"):
                star_info[ANNOTATIONS[3]] += round(
                    float(line.strip("%").split("\t")[1]), 2
                )
            elif line.startswith("% of reads unmapped: too many mismatches"):
                star_info[ANNOTATIONS[4]] = round(
                    float(line.strip("%").split("\t")[1]), 2
                )
            elif line.startswith("% of reads unmapped: too short"):
                star_info[ANNOTATIONS[4]] += round(
                    float(line.strip("%").split("\t")[1]), 2
                )
            elif line.startswith("% of reads unmapped: other"):
                star_info[ANNOTATIONS[4]] += round(
                    float(line.strip("%").split("\t")[1]), 2
                )

    star_info = {key: round(star_info[key], 2) for key in list(star_info.keys())}
    if outfile is None:
        return star_info
    filename = path_leaf(infile)
    filename = filename.strip("Log.final.out")
    counts_df = pd.DataFrame.from_dict(star_info, orient="index").T
    counts_df.index = [filename]
    if outfile:
        counts_df.to_csv(outfile, sep=str("\t"), index=True, header=True)
    return counts_df


def get_strandedness(filepath):
    """Parse output of infer_experiment.py from RSeqC to get strandedness.

    Parameters
    ----------
    filepath : str
               Path to infer_experiment.py output

    Returns
    -------
    strandedness : str
                   reverse or forward or none
    """
    with open(filepath) as f:
        data = f.read()
    splitted = [x.strip() for x in data.split("\n") if len(x.strip()) >= 1]
    assert splitted[0] == "This is SingleEnd Data"
    fwd_percentage = None
    rev_percentage = None
    for line in splitted[1:]:
        if "Fraction of reads failed to determine:" in line:
            continue
        elif 'Fraction of reads explained by "++,--":' in line:
            fwd_percentage = float(line.split(":")[1])
        elif 'Fraction of reads explained by "+-,-+":' in line:
            rev_percentage = float(line.split(":")[1])

    assert rev_percentage is not None
    assert fwd_percentage is not None

    ratio = fwd_percentage / rev_percentage

    if np.isclose([ratio], [1]):
        return "none"
    elif ratio >= 0.5:
        return "forward"
    else:
        return "reverse"


def load_pickle(filepath):
    """Read pickled files easy in Python 2/3"""
    if ".tsv" in filepath:
        raise IndexError
    if sys.version_info > (3, 0):
        pickled = pickle.load(open(filepath, "rb"), encoding="latin1")
    else:
        pickled = pickle.load(open(filepath, "rb"))
    return pickled


def pad_or_truncate(some_list, target_len):
    """Pad or truncate a list upto given target length

    Parameters
    ----------
    some_list : list
                Input list
    target_length : int
                    Final length of list

    If being extended, returns list padded with NAs.
    """
    return some_list[:target_len] + [np.nan] * (target_len - len(some_list))


def pad_five_prime_or_truncate(some_list, offset_5p, target_len):
    """Pad first the 5prime end and then the 3prime end or truncate

    Parameters
    ----------
    some_list : list
                Input list
    offset_5p : int
                5' offset
    target_length : int
                    Final length of list

    If being extended, returns list padded with NAs.
    """
    some_list = list(some_list)
    padded_5p = [np.nan] * offset_5p + some_list
    return padded_5p[:target_len] + [np.nan] * (target_len - len(padded_5p))


def codon_to_anticodon(codon):
    """Codon to anticodon.

    Parameters
    ----------
    codon : string
            Input codon
    """
    pairs = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N"}
    return "".join(pairs[c] for c in codon)[::-1]


def merge_intervals(
    intervals, chromosome_lengths=None, offset_5p=0, offset_3p=0, zero_based=True
):
    """Collapse intervals into non overlapping manner

    Parameters
    ----------
    intervals : list of Interval
    chromosome_lengths : dict
                         A map of each chromosome'e length
                         Only used with offset_3p, offset_5p>0
    offset_5p : int (positive)
                Number of bases to count upstream (5')
    offset_3p : int (positive)
                Number of bases to count downstream (3')
    zero_based: bool
                Indicate if the intervals are zero-based
                True means zero-based half open
                False means one-based full closed

    Returns
    -------
    interval_combined : list of Interval sorted by the start
                        A merged version of intervals
                        This is useful when the annotations are overlapping.
                        Example:
                        chr1 310 320 gene1 +
                        chr1 319 324 gene1 +
                        Returns:
                        chr1 310 324 gene1 +

    gene_offset_5p: Gene wise 5 prime offset
                    This might be different from `offset_5p` in cases where
                    `offset_5p` leads to a negative coordinate
    gene_offset_3p: Gene wise 3 prime offset
                    This might be different from `offset_3p` in cases where
                    `offset_3p` leads to position beyond chromsome length
    """
    if not intervals:
        return ([], offset_5p, offset_3p)

    chroms = list(set([i.chrom for i in intervals]))
    strands = list(set([i.strand for i in intervals]))

    if len(chroms) != 1:
        sys.stderr.write("Error: chromosomes should be unique")
        return ([], offset_5p, offset_3p)
    if len(strands) != 1:
        sys.stderr.write("Error: strands should be unique")
        return ([], offset_5p, offset_3p)

    chrom = chroms[0]
    strand = strands[0]

    # Sort intervals by start
    intervals.sort(key=lambda x: x.start)

    # Find first interval
    first_interval = intervals[0]

    # Find last interval
    last_interval = intervals[-1]
    for i in intervals:
        if i.end > last_interval.end:
            last_interval = i

    if offset_5p != 0 or offset_3p != 0:
        if str(chrom) in chromosome_lengths:
            chrom_length = chromosome_lengths[str(chrom)]
        else:
            warnings.warn("Chromosome {} does not exist".format(chrom), UserWarning)
            chrom_length = np.inf
    else:
        chrom_length = np.inf

    if zero_based:
        lower_bound = 0
    else:
        lower_bound = 1
    upper_bound = chrom_length

    if strand == "+":
        if first_interval.start - offset_5p >= lower_bound:
            first_interval.start -= offset_5p
            gene_offset_5p = offset_5p
        else:
            gene_offset_5p = first_interval.start - lower_bound
            first_interval.start = lower_bound

        if last_interval.end + offset_3p <= upper_bound:
            last_interval.end += offset_3p
            gene_offset_3p = offset_3p
        else:
            gene_offset_3p = upper_bound - last_interval.end
            last_interval.end = upper_bound
    else:
        if last_interval.end + offset_5p <= upper_bound:
            last_interval.end += offset_5p
            gene_offset_5p = offset_5p
        else:
            gene_offset_5p = upper_bound - last_interval.end
            last_interval.end = upper_bound

        if first_interval.start - offset_3p >= lower_bound:
            first_interval.start -= offset_3p
            gene_offset_3p = offset_3p
        else:
            gene_offset_3p = first_interval.start - lower_bound
            first_interval.start = lower_bound

    # Merge overlapping intervals
    to_merge = Interval(chrom, first_interval.start, first_interval.end, strand)
    intervals_combined = []
    for i in intervals:
        if i.start <= to_merge.end:
            to_merge.end = max(to_merge.end, i.end)
        else:
            intervals_combined.append(to_merge)
            to_merge = Interval(chrom, i.start, i.end, strand)
    intervals_combined.append(to_merge)

    return (intervals_combined, gene_offset_5p, gene_offset_3p)


def summarize_counters(samplewise_dict):
    """Summarize gene counts for a collection of samples.

    Parameters
    ----------
    samplewise_dict : dict
                      A dictionary with key as sample name and value
                      as another dictionary of counts for each gene

    Returns
    -------
    totals : dict
             A dictionary with key as sample name and value as total gene count

    """
    totals = {}
    for key, sample_dict in six.iteritems(samplewise_dict):
        totals[key] = np.nansum([np.nansum(d) for d in list(sample_dict.values)])
    return totals


def complementary_strand(strand):
    """Get complementary strand

    Parameters
    ----------
    strand: string
            +/-

    Returns
    -------
    rs: string
        -/+

    """
    if strand == "+":
        return "-"
    elif strand == "-":
        return "+"
    else:
        raise ValueError("Not a valid strand: {}".format(strand))


def read_refseq_bed(filepath):
    """Read refseq bed12 from UCSC.

    Parameters
    ----------
    filepath: string
              Location to bed12

    Returns
    -------
    refseq: dict
            dict with keys as gene name and values as intervaltree

    """
    refseq = defaultdict(IntervalTree)
    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            chrom, tx_start, tx_end, name, score, strand = fields[:6]
            tx_start = int(tx_start)
            tx_end = int(tx_end)
            refseq[chrom].insert(tx_start, tx_end, strand)
    return refseq


def read_bed_as_intervaltree(filepath):
    """Read bed as interval tree

    Useful for reading start/stop codon beds

    Parameters
    ----------
    filepath: string
              Location to bed

    Returns
    -------
    bedint_tree: dict
                 dict with keys as gene name and strand as intervaltree

    """
    bed_df = pybedtools.BedTool(filepath).sort().to_dataframe()
    bed_df["chrom"] = bed_df["chrom"].astype(str)
    bed_df["name"] = bed_df["name"].astype(str)
    bed_grouped = bed_df.groupby("chrom")

    bedint_tree = defaultdict(IntervalTree)
    for chrom, df in bed_grouped:
        df_list = list(zip(df["start"], df["end"], df["strand"]))
        for start, end, strand in df_list:
            bedint_tree[chrom].insert(start, end, strand)
    return bedint_tree


def read_chrom_sizes(filepath):
    """Read chr.sizes file sorted by chromosome name

    Parameters
    ----------
    filepath: string
              Location to chr.sizes

    Returns
    -------
    chrom_lengths: list of tuple
                   A list of tuples with chromsome name and their size
    """

    chrom_lengths = []
    with open(filepath, "r") as fh:
        for line in fh:
            chrom, size = line.strip().split("\t")
            chrom_lengths.append((chrom, int(size)))
        chrom_lengths = list(sorted(chrom_lengths, key=lambda x: x[0]))


def create_bam_index(bam):
    """Create bam index.

    Parameters
    ----------
    bam : str
          Path to bam file
    """
    if isinstance(bam, pysam.AlignmentFile):
        bam = bam.filename
    if not os.path.exists("{}.bai".format(bam)):
        pysam.index(bam)


def is_read_uniq_mapping(read):
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
        nh_count = tags["NH"]
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
            raise RuntimeError("Malformed BAM?")
    if nh_count == 1:
        return True
    return False


def find_first_non_none(positions):
    """Given a list of positions, find the index and value of first non-none element.

    This method is specifically designed for pysam, which has a weird way of returning
    the reference positions. If they are mismatched/softmasked it returns None
    when fetched using get_reference_positions.

    query_alignment_start and query_alignment_end give you indexes of position in the read
    which technically align, but are not softmasked i.e. it is set to None even if the position does not align

    Parameters
    ----------
    positions: list of int
               Positions as returned by pysam.fetch.get_reference_positions

    Return
    ------
    index: int
           Index of first non-None value
    position: int
               Value at that index
    """
    for idx, position in enumerate(positions):
        if position is not None:
            return idx, position


def find_last_non_none(positions):
    """Given a list of positions, find the index and value of last non-none element.


    This function is similar to the `find_first_non_none` function, but does it for the reversed
    list. It is specifically useful for reverse strand cases


    Parameters
    ----------
    positions: list of int
               Positions as returned by pysam.fetch.get_reference_positions

    Return
    ------
    index: int
           Index of first non-None value
    position: int
               Value at that index
    """
    return find_first_non_none(positions[::-1])


# NOTE: We can in principle do a longer metagene anaylsis
# using this helper funciont
def yield_intervals(chrom_size, chunk_size=20000):
    for start in np.arange(0, chrom_size, chunk_size):
        end = start + chunk_size
        if end > chrom_size:
            yield (start, chrom_size)
        else:
            yield (start, end)


def bwsum(bw, chunk_size=5000, scale_to=1e6):
    bw_sum = 0
    if isinstance(bw, six.string_types):
        bw = pyBigWig.open(bw)
    chrom_sizes = bw.chroms()
    for chrom, chrom_size in six.iteritems(chrom_sizes):
        for start, end in yield_intervals(chrom_size, chunk_size):
            bw_sum += np.nansum(bw.values(chrom, start, end))
    scale_factor = 1 / (bw_sum / scale_to)
    return bw_sum, scale_factor


def scale_bigwig(inbigwig, chrom_sizes, outbigwig, scale_factor=1):
    """Scale a bigwig by certain factor.

    Parameters
    ----------
    inbigwig: string
              Path to input bigwig
    chrom_sizes: string
                 Path to chrom.sizes file
    outbigwig: string
               Path to output bigwig
    scale_factor: float
                  Scale by value
    """
    wigfile = os.path.abspath("{}.wig".format(outbigwig))
    chrom_sizes = os.path.abspath(chrom_sizes)
    inbigwig = os.path.abspath(inbigwig)
    outbigwig = os.path.abspath(outbigwig)
    if os.path.isfile(wigfile):
        # wiggletools errors if the file already exists
        os.remove(wigfile)

    cmds = ["wiggletools", "write", wigfile, "scale", str(scale_factor), inbigwig]
    try:
        p = subprocess.Popen(
            cmds,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        stdout, stderr = p.communicate()
        rc = p.returncode
        if rc != 0:
            raise RuntimeError(
                "Error running wiggletools.\nstdout : {} \n stderr : {}".format(
                    stdout, stderr
                )
            )
    except FileNotFoundError:
        raise FileNotFoundError(
            "wiggletool not found on the path." "Use `conda install wiggletools`"
        )

    cmds = ["wigToBigWig", wigfile, chrom_sizes, outbigwig]
    try:
        p = subprocess.Popen(
            cmds,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        stdout, stderr = p.communicate()
        rc = p.returncode
        if rc != 0:
            raise RuntimeError(
                "Error running wigToBigWig.\nstdout : {} \n stderr : {}".format(
                    stdout, stderr
                )
            )
        os.remove(wigfile)
    except FileNotFoundError:
        raise FileNotFoundError(
            "wigToBigwig not found on the path. This is an external "
            "tool from UCSC which can be downloaded from "
            "http://hgdownload.soe.ucsc.edu/admin/exe/. Alternatatively, use "
            "`conda install ucsc-wigtobigwig`"
        )


def get_region_sizes(region_bed):
    """Get summed up size of a CDS/UTR region from bed file

    Parameters
    ----------
    region_bed: string
                Input bed file

    Returns
    -------
    region_sizes: pd.Series
                  Series with region name as index and size as key
    """
    if isinstance(region_bed, six.string_types):
        region_bed = pybedtools.BedTool(region_bed).to_dataframe()
    region_bed_grouped = region_bed.groupby("name")
    region_sizes = {}
    for gene_name, gene_group in region_bed_grouped:
        ## Get rid of trailing dots
        gene_name = re.sub(r"\.[0-9]+", "", gene_name)
        # Collect all intervals at once
        intervals = list(
            zip(
                gene_group["chrom"],
                gene_group["start"],
                gene_group["end"],
                gene_group["strand"],
            )
        )
        for interval in intervals:
            if gene_name not in region_sizes:
                # End is always 1-based so does not require +1
                region_sizes[gene_name] = interval[2] - interval[1]
            else:
                region_sizes[gene_name] += interval[2] - interval[1]
    return pd.Series(region_sizes)


def htseq_to_tpm(htseq_f, outfile, cds_bed_f):
    """Convert htseq-counts file to tpm

    Parameters
    ----------
    htseq_f: string
             Path to htseq-count output
    outfile: string
             Path to output file with tpm values
    cds_bed_f: string
               Path to CDS/genesize bed file
    """
    cds_bed = pybedtools.BedTool(cds_bed_f).to_dataframe()
    cds_bed_sizes = get_region_sizes(cds_bed)
    htseq = pd.read_table(htseq_f, names=["name", "counts"]).set_index("name")
    htseq = htseq.iloc[:-5]
    if htseq.shape[0] <= 10:
        print("Empty dataframe for : {}\n".format(htseq_f))
        return None
    rate = np.log(htseq["counts"]).subtract(np.log(cds_bed_sizes))
    denom = np.log(np.sum(np.exp(rate)))
    tpm = np.exp(rate - denom + np.log(1e6))
    tpm = pd.DataFrame(tpm, columns=["tpm"])
    tpm = tpm.sort_values(by=["tpm"], ascending=False)
    tpm.to_csv(outfile, sep="\t", index=True, header=False)


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


def read_htseq(htseq_f):
    """Read HTSeq file.

    Parameters
    ----------
    htseq_f: str
             Path to htseq counts file

    Returns
    -------
    htseq_df: dataframe
              HTseq counts as in a dataframe
    """
    htseq = pd.read_table(htseq_f, names=["name", "counts"]).set_index("name")
    htseq = htseq.iloc[:-5]
    if htseq.shape[0] <= 10:
        sys.stderr.write("Empty dataframe for : {}\n".format(htseq_f))
        return None
    return htseq


def read_enrichment(
    read_lengths,
    enrichment_range=list(range(28, 33)),
    input_is_stream=False,
    input_is_file=True,
):
    """Calculate read enrichment for a certain range of lengths

    Parameters
    ----------
    read_lengths: Counter
                  A counter with read lengths and their counts
    enrichment_range: range or str
                      Range of reads to concentrate upon
                      (28-32 or range(28,33))
    input_is_stream: bool
                     True if input is sent through stdin

    Returns
    -------
    ratio: float
           Enrichment in this range (Scale 0-1)
    """
    if isinstance(read_lengths, pd.Series):
        pass
    elif input_is_file:
        read_lengths = os.path.abspath(read_lengths)
        if not check_file_exists(read_lengths):
            raise RuntimeError("{} does not exist.".format(read_lengths))
        read_lengths = pd.read_table(read_lengths, sep="\t")
        read_lengths = pd.Series(
            read_lengths["count"].tolist(), index=read_lengths["read_length"].tolist()
        ).add(
            pd.Series([0] * len(enrichment_range), index=list(enrichment_range)),
            fill_value=0,
        )
    elif input_is_stream:
        counter = {}
        for line in read_lengths:
            splitted = list([int(x) for x in line.strip().split("\t")])
            counter[splitted[0]] = splitted[1]
        read_lengths = Counter(counter)
    if isinstance(read_lengths, Counter):
        read_lengths = pd.Series(read_lengths)
        if isinstance(enrichment_range, six.string_types):
            splitted = list([int(x) for x in enrichment_range.strip().split("-")])
        enrichment_range = list(range(splitted[0], splitted[1] + 1))
    rpf_signal = read_lengths.get(list(enrichment_range)).sum()
    total_signal = read_lengths.sum()
    read_lengths = read_lengths.sort_index()
    array = [[x] * int(y) for x, y in zip(read_lengths.index, read_lengths.values)]
    mean_length, std_dev_length = stats.norm.fit(np.concatenate(array).ravel().tolist())

    # mean_length_floor = np.floor(mean_length)
    # 1 - P(x1 < X <x2) = P(X<x1) + P(X>x2) = cdf(x1) + sf(x2)
    cdf_min = stats.norm.cdf(min(enrichment_range), mean_length, std_dev_length)
    sf_max = stats.norm.sf(max(enrichment_range), mean_length, std_dev_length)
    pvalue = cdf_min + sf_max
    ratio = rpf_signal / float(total_signal)
    return ratio, pvalue


def bwshift(bw, shift_by, out_bw, chunk_size=20000):
    """Given a bigwig shift all the values by this
    Shifting by 10:
    variableStep chrom=chr span=1
    1	1
    2	2
    3	5
    4	6
    5	5
    6	3
    7	3
    8	5
    9	5
    10	5
    11	6
    12	6
    13	0
    14	2
    15	3
    16	3
    17	10
    18	4
    19	4
    20	2
    21	2
    22	2
    23	1

    shifted by 10
    variableStep chrom=chr span=1
    1	6
    2	6
    3	0
    4	2
    5	3
    6	3
    7	10
    8	4
    9	4
    10	2
    11	2
    12	2
    13	1
    """
    if isinstance(bw, six.string_types):
        bw = pyBigWig.open(bw)
    chrom_sizes = bw.chroms()
    out_bw = pyBigWig.open(out_bw, "w")
    out_bw.addHeader(list(chrom_sizes.items()), maxZooms=0)
    for chrom, chrom_size in six.iteritems(chrom_sizes):
        for start, end in yield_intervals(chrom_size, chunk_size):
            shifted_values = []
            shifted_end = end + shift_by
            shifted_start = start + shift_by
            # Cases
            # Case1: shifted_start < chrom_start and shifted_end < chrom_end
            # Need to set the first shifted_start - chrom_start bases to 0
            if shifted_start < 0 and shifted_end < chrom_size:
                zero_padding = np.array([0.0] * np.abs(shifted_start))
                shifted_values = np.concatenate(
                    [zero_padding, bw.values(chrom, 0, shifted_end)]
                )
            # Case2: shifted_start < 0 and shifted_end > chrom_end
            elif shifted_start < 0 and shifted_end >= chrom_size:
                zero_padding_start = np.array([0.0] * shifted_start)
                zero_padding_end = np.array([0.0] * (shifted_end - chrom_size))
                shifted_values = np.concatenate(
                    [
                        zero_padding_start,
                        bw.values(chrom, shifted_start, chrom_size),
                        zero_padding_end,
                    ]
                )
            # Case3: shifted_start > 0 and shifted_end > chrom_end
            elif shifted_start > 0 and shifted_end >= chrom_size:
                zero_padding_end = np.array([0.0] * (shifted_end - chrom_size))
                shifted_values = np.concatenate(
                    [bw.values(chrom, shifted_start, chrom_size), zero_padding_end]
                )
            elif shifted_start > 0 and shifted_end < chrom_size:
                shifted_values = bw.values(chrom, shifted_start, shifted_end)
            elif shifted_start >= chrom_size:
                shifted_values = []
            else:
                raise NotImplementedError(
                    "Unhandled case: shifted_start = {} | shifted_end = {} | chrom_size = {}".format(
                        shifted_start, shifted_end, chrom_size
                    )
                )
            assert end - start == len(
                shifted_values
            ), "end = {} | shifted_end = {} | start={} | shifted_start ={} | values_len = {}".format(
                end, shifted_end, start, shifted_start, len(shifted_values)
            )
            starts = np.arange(start, end)
            ends = starts + 1
            out_bw.addEntries([chrom] * (end - start), starts, ends, shifted_values)
