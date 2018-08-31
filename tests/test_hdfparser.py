import pytest
from riboraptor.hdf_parser import HDFParser
from riboraptor.count import export_read_length, export_metagene_coverage
import pandas as pd
import pysam


def test_coverage_from_hdf():
    """Test retrieving coverage.
    Tested by default for 5' orientation.
    """
    region1 = ('chr1', 629938, 634049, '+')
    region2 = ('chr1', 903967, 944088, '+')
    region3 = ('chr1', 903967, 944088, '-')
    fragment_lengths = [24, 25, 21]
    hdf_filepath = 'tests/data/SRX2536403_subsampled.unique.bam_coverage.hdf5'
    hdf = HDFParser(hdf_filepath)
    coverage = hdf.get_coverage(region1, fragment_lengths)
    assert coverage.loc[coverage.start == 629938, '24'].tolist() == [1]
    assert coverage.loc[coverage.start == 634048, '25'].tolist() == [1]

    coverage = hdf.get_coverage(region2, fragment_lengths)
    assert coverage.loc[coverage.start == 944087, '24'].tolist() == [0]
    assert coverage.loc[coverage.start == 944087, '21'].tolist() == [0]

    coverage = hdf.get_coverage(region3, fragment_lengths)
    assert coverage.loc[coverage.start == 944087, '24'].tolist() == [0]
    assert coverage.loc[coverage.start == 944087, '21'].tolist() == [1]
    hdf.close()


def test_length_fragments():
    """Tested if length distribution in hdf
    is same as when obtained from the bam"""
    hdf_filepath = 'tests/data/SRX2536403_subsampled.unique.bam_coverage.hdf5'
    bam_filepath = 'tests/data/SRX2536403_subsampled.unique.bam'
    hdf = HDFParser(hdf_filepath)
    # Ensure everything is int32
    read_lengths_bam = pd.Series(
        export_read_length(bam_filepath)).sort_index().astype('int32')
    read_lengths_hdf = hdf.get_read_length_dist().sort_index().astype('int32')
    assert read_lengths_bam.equals(read_lengths_hdf)


"""

def test_export_coverage():
    bw = 'tests/data/SRX2536403_subsampled.unique.bigWig'
    bed = 'hg38_cds'
    export_metagene_coverage(bed,
                             bw,
                             max_positions=1000,
                             saveto='tests/data/SRX2536403_subsampled.unique.metagene.tsv',
                             offset_5p=0,
                             offset_3p=0)
"""


def test_metagene_coverage():
    #Tested if length distribution in hdf
    #is same as when obtained from the bam
    hdf_filepath = 'tests/data/SRX2536403_subsampled.unique.bam_coverage.hdf5'
    bed = 'hg38_cds'
    hdf = HDFParser(hdf_filepath)
    metagene_coverage_raw, metagene_coverage_normalized = hdf.get_metagene_coverage(
        bed, n_bases=1000)
    assert pd.read_table('tests/data/SRX2536403_subsampled.unique.metagene.tsv'
                         ) == metagene_coverage_normalized['combined_total']


def test_getchromlengths():
    hdf_filepath = 'tests/data/SRX2536403_subsampled.unique.bam_coverage.hdf5'
    bam_filepath = 'tests/data/SRX2536403_subsampled.unique.bam'
    bam = pysam.AlignmentFile(bam_filepath, 'rb')
    reference_and_length = dict(zip(bam.header.references, bam.header.lengths))
    bam.close()
    hdf = HDFParser(hdf_filepath)
    assert hdf.chromosome_lengths == reference_and_length


if __name__ == '__main__':
    pytest.main([__file__])
