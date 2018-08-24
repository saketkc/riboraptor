import h5py
import pandas as pd


class HDFParser(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.h5py_obj = h5py.File(filepath, 'r')

    def close(self):
        self.h5py_obj.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.h5py_obj.close()

    def get_coverage(self,
                     region=(None, None, None, None),
                     fragment_length=None,
                     orientation='5prime'):
        """Get coverage for a selected region.

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

        h5py_fragments = self.h5py_obj['read_lengths']
        if fragment_length == 'all':
            fragment_length = h5py_fragments
        fragment_length = list(map(lambda x: str(x), fragment_length))
        coverages = []
        for l in fragment_length:
            root_obj = self.h5py_obj['fragments'][l][orientation]
            if strand == '-':
                chrom_obj = root_obj['{}_neg'.format(chrom)]
            elif strand == '+':
                chrom_obj = root_obj['{}_pos'.format(chrom)]
            else:
                raise ValueError('strand ill-defined')
            counts_series = pd.Series(
                chrom_obj['counts'], index=chrom_obj['positions'])
            try:
                coverage = counts_series[list(range(start, stop))]
                coverage = coverage.fillna(0).tolist()
            except KeyError:
                # None of the integers in (start, stop) are in
                # the index, so it must be all zero coverage
                coverage = [0] * (stop - start)
            coverages.append(coverage)

        coverages = pd.DataFrame(coverages)
        coverages = coverages.astype(int)
        coverages.columns = list(range(start, stop))
        coverages.index = fragment_length
        coverages = coverages.T.reset_index()
        coverages = coverages.rename(columns={'index': 'start'})
        coverages.columns = ['start'] + list(coverages.columns[1:])
        return coverages

    def get_read_length_dist(self):
        """Get fragment length distribution"""
        return pd.Series(
            self.h5py_obj['read_lengths_counts'],
            index=[int(x) for x in self.h5py_obj['read_lengths']])
