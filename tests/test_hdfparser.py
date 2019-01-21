import pytest
from riboraptor.hdf_parser import HDFParser
from riboraptor.count import export_read_length, export_metagene_coverage
import pandas as pd
import pysam


def test_coverage_from_hdf():
    """Test retrieving coverage.
    Tested by default for 5' orientation.
    """
    region1 = ('chr1', 629916, 629926, '-')
    region2 = ('chr1', 903967, 944088, '+')
    region3 = ('chr1', 903967, 944088, '-')
    fragment_lengths = 'all'  #[24, 25, 21]
    hdf_filepath = 'tests/data/SRX2536403_subsampled.unique.bam_coverage.hdf5'
    hdf = HDFParser(hdf_filepath)
    coverage, coverage_normalized, coverage_sum, coverage_normalized_sum = hdf.get_coverage(
        region1, fragment_lengths, orientation='3prime')
    assert coverage.loc[coverage.start == 629916, '24'].tolist() == [1]
    hdf.close()


def test_length_fragments():
    """Tested if length distribution in hdf
    is same as when obtained from the bam"""
    hdf_filepath = 'tests/data/SRX2536403_subsampled.unique.bam_coverage.hdf5'
    bam_filepath = 'tests/data/SRX2536403_subsampled.unique.bam'
    hdf = HDFParser(hdf_filepath)
    # Ensure everything is int32
    read_lengths_bam = pd.Series(
        export_read_length(bam_filepath)).sort_index().astype('int64')
    read_lengths_hdf = hdf.get_read_length_dist().sort_index().astype('int64')
    assert len(read_lengths_bam) == len(read_lengths_hdf)
    assert list(read_lengths_bam.index) == list(read_lengths_hdf.index)
    assert (read_lengths_bam == read_lengths_hdf).all()


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

"""


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
