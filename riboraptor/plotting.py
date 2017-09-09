"""Plotting methods."""
from itertools import cycle
from itertools import islice

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pandas as pd
import seaborn as sns

from .helpers import round_to_nearest

__FRAME_COLORS__ = ['#1b9e77', '#d95f02', '#7570b3']


def setup_plot():
    """Setup plotting defaults"""
    plt.rcParams['savefig.dpi'] = 120
    plt.rcParams['figure.dpi'] = 120
    plt.rcParams['figure.autolayout'] = False
    plt.rcParams['figure.figsize'] = 12, 8
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['axes.titlesize'] = 20
    plt.rcParams['font.size'] = 10
    plt.rcParams['lines.linewidth'] = 2.0
    plt.rcParams['lines.markersize'] = 8
    plt.rcParams['legend.fontsize'] = 14

    sns.set_style('white')
    sns.set_context('paper', font_scale=2)

def setup_axis(ax, minorticks=10):
    """Setup axes defaults"""
    minor_locator = AutoMinorLocator(minorticks)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.tick_params(which='major', width=2)
    ax.tick_params(which='minor', width=1)
    ax.tick_params(which='major', length=10)
    ax.tick_params(which='minor', length=6)

def plot_framewise_dist(counts, fragment_len_range, ax=None):
    """Plot framewise distribution of fragments.

    Parameters
    ----------
    counts : Series
            A series with position as index and value as counts

    fragment_len_range: int or range
        Range of fragment lengths to average counts over

    ax : matplotlib.Axes
        Default none

    """
    setup_plot()
    assert isinstance(counts, pd.Series)
    if ax is None:
        _, ax = plt.subplots()
    setup_axis(ax)
    ax.set_ylabel('Number of reads')
    ax.set_xlim(min(counts.index) - 0.6,
                round_to_nearest(max(counts.index), 10) + 0.6)
    barlist = ax.bar(counts.index, counts.values)
    barplot_colors = list(islice(cycle(__FRAME_COLORS__), None, len(counts.index)))
    for index, cbar in enumerate(barlist):
        cbar.set_color(barplot_colors[index])
    ax.legend((barlist[0], barlist[1], barlist[2]),
              ('Frame 1', 'Frame 2', 'Frame 3'))

    sns.despine(trim=True, offset=20)
    return ax


def plot_fragment_dist(fragment_lengths, ax=None):
    """Plot fragment length distribution.

    Parameters
    ----------
    fragment_lengths : array_like
                     Array of fragment lengths

    ax : matplotlib.Axes
        Axis object

    """
    setup_plot()
    if ax is None:
        _, ax = plt.subplots()
    setup_axis(ax, 5)
    fragment_lengths = pd.Series(fragment_lengths)
    fragment_lengths_counts = fragment_lengths.value_counts().sort_index()

    ax.bar(fragment_lengths_counts.index, fragment_lengths_counts)
    ax.set_xlim(min(fragment_lengths_counts.index) - 0.5,
                round_to_nearest(max(fragment_lengths_counts.index), 10) + 0.5)
    sns.despine(trim=True, offset=20)
    return ax

def plot_rpf_density(counts, ax=None, marker=False, color='royalblue', label=None):
    """Plot RPF density around start/stop codons.

    Parameters
    ----------
    counts : Series
        A series with coordinates as index and counts as values
    ax : matplotlib.Axes
        Axis to create object on
    marker : string
        'o'/'x'
    color : string
        Line color
    label : string
        Label (useful only if plotting multiple objects on same axes)

    """
    setup_plot()
    if ax is None:
        _, ax = plt.subplots()
    setup_axis(ax)
    ax.set_ylabel('Number of reads')
    if not marker:
        ax.plot(counts.index, counts.values, color=color,
                linewidth=2, label=label)
    else:
        ax.plot(counts.index, counts.values, color=color,
                marker='o', linewidth=2, label=label)
    ax.set_xlim(round_to_nearest(ax.get_xlim()[0], 50) - 0.6,
                round_to_nearest(ax.get_xlim()[1], 50) + 0.6)
    sns.despine(trim=True, offset=20)
    return ax
