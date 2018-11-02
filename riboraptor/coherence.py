from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import pickle
import numpy as np
import pandas as pd
import six
from scipy import signal
from .helpers import identify_peaks
from .statistics import coherence_pvalue


def _shift_bit_length(x):
    """Shift bit"""
    return 1 << (x - 1).bit_length()


def _padwithzeros(vector, pad_width, iaxis, kwargs):
    """Pad with zeros"""
    vector[:pad_width[0]] = 0
    vector[-pad_width[1]:] = 0
    return vector


def naive_periodicity(values, identify_peak=False):
    '''Calculate periodicity in a naive manner

    Take ratio of frame1 over avg(frame2+frame3) counts. By default
    the first value is treated as the first frame as well

    Parameters
    ----------
    values : Series
             Metagene profile

    Returns
    -------
    periodicity : float
                  Periodicity

    '''
    if identify_peak:
        peak_location = identify_peaks(values)
        min_index = min(values.index)
        max_index = max(values.index)
        # Frame 1 starts at the peak_location
        # We assume that if values is series it will have continuous
        # indices
        values = values[np.arange(peak_location, max_index)]
    values = pd.Series(values)
    values = values.fillna(0)
    frame1_total = 0
    frame2_total = 0
    frame3_total = 0
    values = list(values)[0:len(values) - len(values) % 3]
    for i in range(0, len(values), 3):
        frame1_total += values[i]
    for i in range(1, len(values), 3):
        frame2_total += values[i]
    for i in range(2, len(values), 3):
        frame3_total += values[i]
    return frame1_total / (frame1_total + frame2_total + frame3_total)


def coherence_ft(values, nperseg=30, noverlap=15, window='flattop'):
    """Calculate coherence and an idea ribo-seq signal based on Fourier transform

    Parameters
    ----------
    values : array like
             List of values

    Returns
    -------
    periodicity : float
                  Periodicity score calculated as
                  coherence between input and idea 1-0-0 signal

    window: str or tuple
            See scipy.signal.get_window


    f: array like
       List of frequencies

    Cxy: array like
         List of coherence at the above frequencies

    """
    length = len(values)
    uniform_signal = [1, 0, 0] * (length // 3)
    mean_centered_values = values - np.nanmean(values)
    normalized_values = mean_centered_values / \
        np.max(np.abs(mean_centered_values))

    mean_centered_values = uniform_signal - np.nanmean(uniform_signal)
    uniform_signal = mean_centered_values / \
        np.max(np.abs(uniform_signal))
    f, Cxy = signal.coherence(
        normalized_values,
        uniform_signal,
        window=window,
        nperseg=nperseg,
        noverlap=noverlap)
    periodicity_score = Cxy[np.argwhere(np.isclose(f, 1 / 3.0))[0]][0]
    return periodicity_score, f, Cxy


def coherence(original_values):
    """Calculate coherence and an idea ribo-seq signal

    Parameters
    ----------
    values : array like
             List of values

    Returns
    -------
    periodicity : float
                  Periodicity score calculated as
                  coherence between input and idea 1-0-0 signal

    f: array like
       List of frequencies

    Cxy: array like
         List of coherence at the above frequencies

    """
    if not isinstance(original_values, list):
        original_values = list(original_values)
    coh, pval, valid = 0.0, 1.0, -1
    for frame in [0, 1, 2]:
        values = original_values[frame:]
        normalized_values = []
        i = 0
        while i + 2 < len(values):
            if values[i] == values[i + 1] == values[i + 2] == 0:
                i += 3
                continue
            real = values[i] + values[i + 1] * np.cos(
                2 * np.pi / 3) + values[i + 2] * np.cos(4 * np.pi / 3)
            imag = values[i + 1] * np.sin(
                2 * np.pi / 3) + values[i + 2] * np.sin(4 * np.pi / 3)
            norm = np.sqrt(real**2 + imag**2)
            if norm == 0:
                norm = 1
            normalized_values += [
                values[i] / norm, values[i + 1] / norm, values[i + 2] / norm
            ]
            i += 3

        length = len(normalized_values) // 3 * 3
        if length == 0:
            return (0.0, 1.0, 0)
        normalized_values = normalized_values[:length]
        uniform_signal = [1, 0, 0] * (len(normalized_values) // 3)
        f, Cxy = signal.coherence(
            normalized_values,
            uniform_signal,
            window=[1.0, 1.0, 1.0],
            nperseg=3,
            noverlap=0)
        try:
            periodicity_score = Cxy[np.argwhere(np.isclose(f, 1 / 3.0))[0]][0]
            periodicity_pval = coherence_pvalue(periodicity_score, length // 3)
        except:
            periodicity_score = 0.0
            periodicity_pval = 1.0
        if periodicity_score > coh:
            coh = periodicity_score
            pval = periodicity_pval
            valid = length
        if valid == -1:
            valid = length
    return coh, pval, valid


def get_periodicity(values, input_is_stream=False):
    """Calculate periodicty wrt 1-0-0 signal.

    Parameters
    ----------
    values : array like
             List of values

    Returns
    -------
    periodicity : float
                  Periodicity calculated as cross
                  correlation between input and idea 1-0-0 signal
    """
    from mtspec import mtspec, mt_coherence
    tbp = 4
    kspec = 3
    nf = 30
    p = 0.90
    if input_is_stream:
        values = list(map(lambda x: float(x.rstrip()), values))
    if isinstance(values, six.string_types):
        try:
            values = pickle.load(open(values))
        except KeyError:
            pass
    values = pd.Series(values)
    values = values[0:max(values.index)]
    length = len(values)
    next_pow2_length = _shift_bit_length(length)
    values = np.lib.pad(values,
                        (0, next_pow2_length - len(values) % next_pow2_length),
                        _padwithzeros)
    mean_centered_values = values - np.nanmean(values)
    normalized_values = mean_centered_values / \
        np.max(np.abs(mean_centered_values))
    uniform_signal = [1, -0.6, -0.4] * (next_pow2_length // 3)
    uniform_signal = np.lib.pad(
        uniform_signal,
        (0, next_pow2_length - len(uniform_signal) % next_pow2_length),
        _padwithzeros)
    out = mt_coherence(
        1,
        normalized_values,
        uniform_signal,
        tbp,
        kspec,
        nf,
        p,
        freq=True,
        phase=True,
        cohe=True,
        iadapt=1)
    spec, freq, jackknife, fstatistics, _ = mtspec(
        data=values,
        delta=1.,
        time_bandwidth=4,
        number_of_tapers=3,
        statistics=True,
        rshape=0,
        fcrit=0.9)
    p = 99
    base_fstat = np.percentile(fstatistics, p)
    fstat = fstatistics[np.argmin(np.abs(freq - 1 / 3.0))]
    p_value = 'p < 0.05'
    if fstat < base_fstat:
        p_value = 'p > 0.05'
    return out['cohe'][np.argmin(np.abs(out['freq'] - 1 / 3.0))], p_value
