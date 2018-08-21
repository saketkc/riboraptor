import h5py
import pandas as pd


class HDFParser(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.h5py_obj = h5py.File(filepath, 'r')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.h5py_obj.close()

    def get_coverage(self, region=(None, None, None, None), fragment_length=None):
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
        assert chrom is not None, 'chromosome not set'
        assert start is not None, 'start position is not set'
        assert stop is not None, 'end position not set'
        assert strand is not None, 'strand is not set'
        assert fragment_length is not None, 'fragment_length not set'

        if isinstance(fragment_length, int):
            fragment_length = [fragment_length]

        h5py_fragments = self.h5py_obj.keys()
        if fragment_length == 'all':
            fragment_length = h5py_fragments

        coverages = []
        for l in fragment_length:
            root_obj = self.h5py_obj[l]
            if strand == '-':
                chrom_obj = root_obj['{}_neg'.format(chrom)]
            elif strand == '+':
                chrom_obj = root_obj['{}_pos'.format(chrom)]
            else:
                raise ValueError('strand ill-defined')
            #chrom_size = chrom_obj['chrom_size'][0]
            #null_series = pd.Series([0]*chrom_size, index=range(chrom_size))
            counts_series = pd.Series(chrom_obj['counts'], index=chrom_obj['positions'])
            try:
                coverage = counts_series[range(start, stop)]
                coverage = coverage.fillna(0).tolist()
            except KeyError:
                # None of the integers in (start, stop) are in
                # the index, so it must be all zero coverage
                coverage =  [0]*(stop-start)
            coverages.append(coverage)

        coverages = pd.DataFrame(coverages, header=None)
        coverages.columns = range(start, stop)
        coverages.index = fragment_length
        return coverages




