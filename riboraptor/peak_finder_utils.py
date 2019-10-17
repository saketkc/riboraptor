import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy import signal
from scipy import stats
from scipy.signal import savgol_filter


def collpase_gene_coverage_to_codon(gene_profile):
    """Assume the gene 0 is the true zero and exclude
    trailing positions which are not 0 mod 3
    """
    codon_profile = []
    gene_profile = np.array(gene_profile)
    for i in range(0, len(gene_profile) - 3, 3):
        codon_profile.append(np.nansum(gene_profile[np.arange(i, i + 3)]))
    # return pd.Series(codon_profile, index=np.arange(1, len(codon_profile)+1))
    return codon_profile  # , index=np.arange(1, len(codon_profile)+1))


def g_transform(data):
    data = np.array(data)
    return np.log(np.log(data + 1) + 1)


def inverse_g_transform(y):
    y = np.array(y)
    return np.exp(np.exp(y) - 1) - 1


def baseline_correct(y, n_iterations=100):
    z = g_transform(y)
    z_copy = np.empty(len(z))
    n = len(z)
    for i in range(0, n_iterations):
        for j in range(i, n - i):
            mean_z = 0.5 * (z[j - i] + [j + i])
            mean_z = min(mean_z, z[j])
            z_copy[j] = mean_z
        for k in range(i, n - i):
            z[k] = z_copy[k]
    inv_z = inverse_g_transform(z)
    return inv_z


def med_abs_dev(data):
    """Calculate Median absolute deviation
    """
    return 1.4826 * max(np.nanmedian(np.abs(data - np.nanmedian(data))), 1e-4)


def calculate_peaks(data, order=3, snr=2.5):
    """ Calculate Peaks
    """
    if isinstance(data, pd.Series):
        index = data.index
    else:
        index = np.arange(0, len(data))
    data = np.array(data)
    data_rel_max_idx = signal.argrelmax(data, axis=0, order=order)[0]
    noise = med_abs_dev(data)
    # peaks_height = np.zeros(len(data))
    peaks_idx = [x for x in data_rel_max_idx if data[x] > snr * noise]
    peaks_x = index[peaks_idx]
    # for x in peaks_idx:
    #    peaks_height[x] = data[x]
    peaks_height = data[peaks_idx]
    return list(peaks_x), list(peaks_height)


def calculate_snr(data):
    data = np.array(data)
    sigma = med_abs_dev(data)
    return data / sigma


def gaussian_pvalue(values):
    values = np.array(values)
    mean = np.nanmean(values)
    std = np.std(values)
    zscore = (values - mean) / std
    pvalues = stats.norm.sf(zscore)
    return pvalues


def plot_orf_with_peak(profile, peaks_location):
    peaks_height = np.array(profile)[peaks_location]

    filtered = savgol_filter(profile, 15, 3)
    fig = plt.figure(figsize=(8, 8))
    ax1 = plt.subplot(311)
    ax1.plot(profile)
    ax1.set_ylabel("RPF count")

    ax2 = plt.subplot(312, sharex=ax1)
    peaks_height_exp = pd.Series([0] * len((profile)))
    peaks_height_exp[peaks_location] = peaks_height
    ax2.stem(peaks_height_exp.index, peaks_height_exp)
    ax2.set_ylabel("RPF count")

    ax3 = plt.subplot(313, sharex=ax1)
    pvalues = -np.log10(gaussian_pvalue(filtered))
    ax3.plot(pvalues)
    ax3.axhline(y=2, color="#D55E00", linestyle="dashed")
    ax3.set_ylabel("-log(pval)")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    return fig
