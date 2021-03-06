"""Command line interface for riboraptor
"""

import os
import sys
from textwrap import dedent

import click
import glob
from click_help_colors import HelpColorsGroup
import six
import pandas as pd

from . import __version__

from .count import export_gene_coverages
from .count import export_metagene_coverage
from .count import export_read_counts
from .count import merge_gene_coverages
from .count import merge_read_counts
from .count import export_read_length
from .count import read_enrichment
from .count import bedgraph_to_bigwig
from .count import bam_to_bedgraph
from .count import count_uniq_mapping_reads
from .count import extract_uniq_mapping_reads
from .count import get_bam_coverage
from .count import get_bam_coverage_on_bed

from .helpers import path_leaf
from .helpers import scale_bigwig
from .loomify import write_loom_batches
from .infer_protocol import infer_protocol

from .sequence import export_gene_sequences

from .download import run_download_sra_script
from .coherence import naive_periodicity
from .plotting import plot_read_counts
from .plotting import plot_read_length_dist
from .plotting import plot_periodicity_df
from .hdf_parser import create_metagene_from_multi_bigwig
from .hdf_parser import hdf_to_bigwig
from .hdf_parser import tsv_to_bigwig
from .hdf_parser import merge_bigwigs
from .hdf_parser import HDFParser
from .hdf_parser import normalize_bw_hdf
from .helpers import bwsum
from .helpers import bwshift
from .orf_seq import orf_seq

click.disable_unicode_literals_warning = True
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(
    cls=HelpColorsGroup, help_headers_color="yellow", help_options_color="green"
)
@click.version_option(version=__version__)
def cli():
    """riboraptor: Tool for ribosome profiling analysis"""
    pass


###################### export-gene-coverages #################################
@cli.command(
    "export-gene-coverages",
    context_settings=CONTEXT_SETTINGS,
    help="Export gene level coverage for all genes of given region",
)
@click.option(
    "--bed",
    help="Path to bed file or a genome name (hg38_utr5, hg38_cds, hg38_utr3)",
    required=True,
)
@click.option("--bw", help="Path to bigwig file", required=True)
@click.option("--saveto", help="Path to write output", default=None, show_default=True)
@click.option(
    "--offset_5p",
    help="Number of upstream bases to count(5')",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "--offset_3p",
    help="Number of downstream bases to count(3')",
    type=int,
    default=0,
    show_default=True,
)
def export_gene_coverages_cmd(bed, bw, saveto, offset_5p, offset_3p):
    export_gene_coverages(bed, bw, saveto, offset_5p, offset_3p)


###################### export-metagene-coverages ##############################
@cli.command(
    "export-metagene-coverage",
    context_settings=CONTEXT_SETTINGS,
    help="Export metagene coverage for given region",
)
@click.option(
    "--bed",
    help="Path to bed file or a genome name (hg38_utr5, hg38_cds, hg38_utr3)",
    required=True,
)
@click.option("--bw", help="Path to bigwig file", required=True)
@click.option(
    "--max_positions",
    help="maximum positions to count",
    type=int,
    default=500,
    show_default=True,
)
@click.option("--saveto", help="Path to write output", default=None, show_default=True)
@click.option(
    "--offset_5p",
    help="Number of upstream bases to count(5')",
    type=int,
    default=60,
    show_default=True,
)
@click.option(
    "--offset_3p",
    help="Number of downstream bases to count(3')",
    type=int,
    default=60,
    show_default=True,
)
@click.option(
    "--orientation",
    help="Orientation of read ends",
    type=click.Choice(["5prime", "3prime"]),
    required=True,
)
@click.option(
    "--threads",
    help="Number of threads to use",
    type=int,
    default=16,
    show_default=True,
)
def export_metagene_coverage_cmd(
    bed, bw, max_positions, saveto, offset_5p, offset_3p, orientation, threads
):
    metagene_profile = export_metagene_coverage(
        bed, bw, max_positions, saveto, offset_5p, offset_3p, orientation, threads
    )
    if saveto is None:
        for i, count in six.iteritems(metagene_profile):
            sys.stdout.write("{}\t{}".format(i, count))
            sys.stdout.write(os.linesep)


###################### export-read-counts ##############################
@cli.command(
    "export-read-counts",
    context_settings=CONTEXT_SETTINGS,
    help="Export read counts from gene coverages file",
)
@click.option("--gene_coverages", help="Path to gene coverages file", required=True)
@click.option("--saveto", help="Path to write output", default=None, show_default=True)
@click.option("--keep_offsets", help="whether keep the 5' and 3' offsets", is_flag=True)
def export_read_counts_cmd(gene_coverages, saveto, keep_offsets):
    export_read_counts(gene_coverages, saveto, keep_offsets)


###################### merge-gene-coverages ##############################
@cli.command(
    "merge-gene-coverages",
    context_settings=CONTEXT_SETTINGS,
    help="merge gene coverages to generate metagene coverage",
)
@click.option("--gene_coverages", help="Path to gene coverages file", required=True)
@click.option(
    "--max_positions",
    help="maximum positions to count",
    type=int,
    default=500,
    show_default=True,
)
@click.option("--saveto", help="Path to write output", default=None, show_default=True)
def merge_gene_coverages_cmd(gene_coverages, max_positions, saveto):
    merge_gene_coverages(gene_coverages, max_positions, saveto)


###################### merge-read-counts ##############################
@cli.command(
    "merge-read-counts",
    context_settings=CONTEXT_SETTINGS,
    help="merge read counts to generate count table",
)
@click.option(
    "--read_counts", help="Path to file containing read counts paths", required=True
)
@click.option("--saveto", help="Path to write output", required=True)
def merge_read_counts_cmd(read_counts, saveto):
    merge_read_counts(read_counts, saveto)


#################### export-read-length ######################################
@cli.command(
    "export-read-length",
    context_settings=CONTEXT_SETTINGS,
    help="Calculate read length distribution",
)
@click.option("--bam", help="Path to BAM file", required=True)
@click.option("--saveto", help="Path to write read length dist tsv output")
def export_read_length_cmd(bam, saveto):
    counts = export_read_length(bam, saveto)
    for i, count in six.iteritems(dict(counts)):
        sys.stdout.write("{}\t{}".format(i, count))
        sys.stdout.write(os.linesep)


#################### read-enrichment ######################################
@cli.command(
    "read-enrichment",
    context_settings=CONTEXT_SETTINGS,
    help="Calculate read enrichment for a certain range of lengths",
)
@click.option("--read_lengths", help="Path to read length tsv", required=True)
@click.option("--min_length", help="The low end of the range", default=28)
@click.option("--max_length", help="The high end of the range", default=32)
def read_enrichment_cmd(read_lengths, min_length, max_length):
    ratio = read_enrichment(read_lengths, min_length, max_length)
    sys.stdout.write(
        "Enrichment of length range {}-{}: {}".format(min_length, max_length, ratio)
    )
    sys.stdout.write(os.linesep)


#################### periodicity #############################################
@cli.command(
    "periodicity", context_settings=CONTEXT_SETTINGS, help="Calculate periodicity"
)
@click.option("--counts", help="Path to counts file (if not stdin)")
def periodicity_cmd(counts):
    if counts:
        counts = pd.read_table(counts)
        counts = pd.Series(counts["count"].tolist(), index=counts["position"].tolist())
        periodicity = naive_periodicity(counts)
    else:
        periodicity = naive_periodicity(sys.stdin.readlines(), input_is_stream=True)
    sys.stdout.write("Periodicity: {}".format(periodicity))
    sys.stdout.write(os.linesep)


#################### export-gene-sequences ####################################
@cli.command(
    "export-gene-sequences",
    context_settings=CONTEXT_SETTINGS,
    help="Export gene level sequence for all genes of given region",
)
@click.option(
    "--bed",
    help="Path to bed file or a genome name (hg38_utr5, hg38_cds, hg38_utr3)",
    required=True,
)
@click.option("--fasta", help="Path to fasta file", required=True)
@click.option("--saveto", help="Path to write output", default=None, show_default=True)
@click.option(
    "--offset_5p", help="Number of upstream bases to count(5')", type=int, default=0
)
@click.option(
    "--offset_3p", help="Number of downstream bases to count(3')", type=int, default=0
)
def export_gene_sequences_cmd(bed, fasta, saveto, offset_5p, offset_3p):
    export_gene_sequences(bed, fasta, saveto, offset_5p, offset_3p)


######################## plot-metagene ####################################
@cli.command(
    "plot-metagene", context_settings=CONTEXT_SETTINGS, help="Plot metagene profile"
)
@click.option("--counts", help="Path to counts file (if not stdin)")
@click.option("--title", help="Plot Title")
@click.option(
    "--marker",
    help="Marker (o/x) for plots",
    type=click.Choice(["o", "x"]),
    default=None,
)
@click.option("--color", help="Color", default="royalblue")
@click.option(
    "--millify_labels", help="Convert labels on Y-axis to concise form?", is_flag=True
)
@click.option("--identify_peak", help="Identify Peak?", is_flag=True)
@click.option("--positions", help="Range of positions to plot", default="-60:500")
@click.option(
    "--saveto",
    help="Path to file (png/pdf) to save to",
    default=None,
    show_default=True,
)
@click.option("--ylabel", help="Y axix label", default="Normalized RPF density")
@click.option("--yrotation", default=0, help="Rotate y axis labels by", type=int)
@click.option("--xrotation", default=0, help="Rotate x axis labels by", type=int)
@click.option("--ascii", help="Plot ascii", is_flag=True)
def plot_read_counts_cmd(
    counts,
    title,
    marker,
    color,
    millify_labels,
    identify_peak,
    positions,
    saveto,
    ylabel,
    xrotation,
    yrotation,
    ascii,
):
    if counts:
        plot_read_counts(
            counts,
            title=title,
            marker=marker,
            color=color,
            millify_labels=millify_labels,
            position_range=positions,
            identify_peak=identify_peak,
            ylabel=ylabel,
            saveto=saveto,
            ascii=ascii,
            xrotation=xrotation,
            yrotation=yrotation,
        )
    else:
        plot_read_counts(
            sys.stdin.readlines(),
            title=title,
            marker=marker,
            color=color,
            millify_labels=millify_labels,
            identify_peak=identify_peak,
            position_range=positions,
            saveto=saveto,
            ylabel=ylabel,
            ascii=ascii,
            input_is_stream=True,
            xrotation=xrotation,
            yrotation=yrotation,
        )


######################## plot-read-length ####################################
@cli.command(
    "plot-read-length",
    context_settings=CONTEXT_SETTINGS,
    help="Plot read length distribution",
)
@click.option("--read-lengths", help="Path to read length counts file")
@click.option("--hdf", help="Path to hdf file")
@click.option("--title", help="Plot Title")
@click.option(
    "--millify_labels", help="Convert labels on Y-axis to concise form?", is_flag=True
)
@click.option(
    "--saveto",
    help="Path to file (png/pdf) to save to",
    default=None,
    show_default=True,
)
@click.option("--ascii", help="Do not plot ascii", is_flag=True)
def plot_read_length_dist_cmd(read_lengths, hdf, title, millify_labels, saveto, ascii):
    if read_lengths:
        plot_read_length_dist(
            read_lengths,
            title=title,
            millify_labels=millify_labels,
            input_is_stream=False,
            saveto=saveto,
            ascii=ascii,
        )
    elif hdf:
        path, ext = os.path.splitext(saveto)
        h5 = HDFParser(hdf)
        read_length_dist = h5.get_read_length_dist()
        query_alignment_length_dist = h5.get_query_alignment_length_dist()
        plot_read_length_dist(
            read_length_dist,
            title=title,
            millify_labels=millify_labels,
            input_is_stream=False,
            saveto=saveto,
            ascii=ascii,
        )
        plot_read_length_dist(
            query_alignment_length_dist,
            title=title,
            millify_labels=millify_labels,
            input_is_stream=False,
            saveto="{}_query_aligned_lengths{}".format(path, ext),
            ascii=ascii,
        )
        h5.close()
    else:
        plot_read_length_dist(
            sys.stdin.readlines(),
            title=title,
            millify_labels=millify_labels,
            input_is_stream=True,
            saveto=saveto,
            ascii=ascii,
        )


####################### bam-to-bedgraph ######################################
@cli.command(
    "bam-to-bedgraph", context_settings=CONTEXT_SETTINGS, help="Convert bam to bedgraph"
)
@click.option("--bam", "-i", help="Path to BAM file", required=True)
@click.option(
    "--strand",
    "-s",
    help="Count from strand of this type only",
    type=click.Choice(["+", "-", "both"]),
    default="both",
)
@click.option(
    "--end_type",
    "-e",
    help="Pileup 5' / 3'/ either ends",
    type=click.Choice(["5prime", "3prime", "either"]),
    default="5prime",
)
@click.option(
    "--saveto",
    "-o",
    help="Path to write bedgraph output",
    default=None,
    show_default=True,
)
def bam_to_bedgraph_cmd(bam, strand, end_type, saveto):
    bedgraph = bam_to_bedgraph(bam, strand, end_type, saveto)
    if saveto is None:
        sys.stdout.write(str(bedgraph))
        sys.stdout.write(os.linesep)


####################### uniq-bam ##############################################
@cli.command(
    "uniq-bam",
    context_settings=CONTEXT_SETTINGS,
    help="Create a new bam with unique mapping reads only",
)
@click.option("--inbam", required=True)
@click.option("--outbam", required=True)
def extract_uniq_mapping_reads_cmd(inbam, outbam):
    extract_uniq_mapping_reads(inbam, outbam)


####################### bedgraph-to-bigwig ######################################
@cli.command(
    "bedgraph-to-bigwig",
    context_settings=CONTEXT_SETTINGS,
    help="Convert bedgraph to bigwig",
)
@click.option(
    "--bedgraph",
    "-bg",
    "-i",
    help="Path to bedgraph file (optional)",
    default=None,
    show_default=True,
)
@click.option("--sizes", "-s", help="Path to genome chrom.sizes file", required=True)
@click.option("--saveto", "-o", help="Path to write bigwig output", required=True)
def bedgraph_to_bigwig_cmd(bedgraph, sizes, saveto):
    if bedgraph:
        bedgraph_to_bigwig(bedgraph, sizes, saveto)
    else:
        bedgraph_to_bigwig(sys.stdin.readlines(), sizes, saveto, True)


###################### uniq-mapping-count ######################################
@cli.command(
    "uniq-mapping-count",
    context_settings=CONTEXT_SETTINGS,
    help="Count number of unique mapping reads",
)
@click.option("--bam", help="Path to BAM file", required=True)
def uniq_mapping_cmd(bam):
    count = count_uniq_mapping_reads(bam)
    sys.stdout.write(str(count))
    sys.stdout.write(os.linesep)


###################### get-bam-metagene-coverage ######################################
@cli.command(
    "bam-metagene-coverage",
    context_settings=CONTEXT_SETTINGS,
    help="Get metagene coverage from bam",
)
@click.option("--bam", help="Path to BAM file", required=True)
@click.option("--bed", help="Path to bed file", required=True)
@click.option(
    "--protocol", default="stranded", help="Library preparation protocol", required=True
)
@click.option("--orientation", default="5prime", help="track 5prime or 3prime")
@click.option(
    "--max_positions",
    help="maximum positions to count",
    type=int,
    default=1000,
    show_default=True,
)
@click.option(
    "--offset", help="Number of upstream bases to count", type=int, default=60
)
@click.option("--saveto", help="Path to store coverage stats", required=True)
def get_bam_coverage_on_bed_cmd(
    bam, bed, protocol, orientation, max_positions, offset, saveto
):
    get_bam_coverage_on_bed(
        bam, bed, protocol, orientation, max_positions, offset, saveto
    )


###################### download function #########################################
@cli.command(
    "download-sra", context_settings=CONTEXT_SETTINGS, help="Download SRA data"
)
@click.option("--out", help="root directory to download all datasets")
@click.option("-a", "--ascp", help="Path to ascp private key")
@click.option(
    "-f", "--file", "--srpfile", help="File containing list of SRPs one per line"
)
@click.argument("srp_id_list", nargs=-1, required=True)
def download_srp_cmd(out, ascp, srpfile, srp_id_list):
    run_download_sra_script(out, ascp, srpfile, list(srp_id_list))


###################### infer-protocol function #########################################
@cli.command(
    "infer-protocol", context_settings=CONTEXT_SETTINGS, help="Infer protocol from BAM"
)
@click.option("--bam", help="Path to bam file")
@click.option(
    "--refseq", help="Path to reseq file to be used for defining the gene strands"
)
@click.option(
    "--n_reads",
    type=int,
    default=10000,
    help="Number of mapped reads to use for estimation",
)
def infer_protocol_cmd(bam, refseq, n_reads):
    protocol, forward_mapped, reverse_mapped, total = infer_protocol(
        bam, refseq, n_reads
    )
    print(
        dedent(
            """\
                 Forward mapped proportion: {:.4f}
                 Reverse mapped proportion: {:.4f}
                 Likely protocol: {}""".format(
                forward_mapped, reverse_mapped, protocol
            )
        )
    )


################### Create metagene from multiple bigwigs #####################
@cli.command(
    "metagene-multi-bw",
    context_settings=CONTEXT_SETTINGS,
    help="Merge multiple bigwigs",
)
@click.option("--pattern", type=str, help="Pattern to search for", required=True)
@click.option("--bed", type=str, help="Path to BED file", required=True)
@click.option("--saveto", type=str, help="Save output bigwig to", required=True)
@click.option(
    "--max_positions",
    help="maximum positions to count",
    type=int,
    default=500,
    show_default=True,
)
@click.option(
    "--offset_5p",
    help="Number of upstream bases to count(5')",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "--offset_3p",
    help="Number of downstream bases to count(3')",
    type=int,
    default=0,
    show_default=True,
)
def merge_multiple_bw(pattern, bed, saveto, max_positions, offset_5p, offset_3p):
    bigwigs = glob.glob(pattern, recursive=True)
    create_metagene_from_multi_bigwig(
        bed, bigwigs, max_positions, offset_5p, offset_3p, n_jobs=16, saveto=saveto
    )


##################### Create bigwig from tsv #################################
@cli.command(
    "tsv-to-bw", context_settings=CONTEXT_SETTINGS, help="Create bigwig from tsv"
)
@click.option("--tsv", type=str, help="Path to tsv file", required=True)
@click.option("--chromsizes", type=str, help="Path to chrom.sizes file", required=True)
@click.option("--prefix", type=str, help="Prefix ", required=True)
def tsv_to_bw_cmd(tsv, chromsizes, prefix):
    tsv_to_bigwig(tsv, chromsizes, prefix)


################### Create metagene from multiple bigwigs #####################
@cli.command(
    "hdf-to-bw", context_settings=CONTEXT_SETTINGS, help="Create bigwig from hdf"
)
@click.option("--hdf", type=str, help="Path to hdf file", required=True)
@click.option("--prefix", type=str, help="Prefix ", required=True)
@click.option(
    "--readlength",
    type=int,
    help="Create bw only of this fragment length",
    required=False,
    default="all",
)
def hdf_to_bw_cmd(hdf, prefix, readlength):
    hdf_to_bigwig(hdf, prefix, readlength)


###################### get-bam-coverage ######################################
@cli.command(
    "bam-coverage",
    context_settings=CONTEXT_SETTINGS,
    help="Get strandwise coerage from bam",
)
@click.option("--bam", help="Path to BAM file", required=True)
@click.option("--genebed", help="Path to genes.bed file", required=True)
@click.option("--outprefix", help="Prefix to store coverage output", required=True)
def bam_coverage_cmd(bam, genebed, outprefix):
    get_bam_coverage(bam, genebed, outprefix)


################### Merge multiple bigwigs ####################################
@cli.command(
    "merge-bw", context_settings=CONTEXT_SETTINGS, help="Merge multiple bigwigs"
)
@click.option("--pattern", type=str, help="Pattern to search for", required=True)
@click.option("--chromsizes", type=str, help="Path to chrom.sizes file", required=True)
@click.option("--saveto", type=str, help="Save output bigwig to", required=True)
def merge_bw_cmd(pattern, chromsizes, saveto):
    bigwigs = glob.glob(os.path.abspath(pattern), recursive=True)
    merge_bigwigs(bigwigs, chromsizes, saveto)


####################### Sum bigwig ##########################################
@cli.command("bwsum", context_settings=CONTEXT_SETTINGS, help="Scale bigwig")
@click.option("--inbw", type=str, help="Path to input bigwig", required=True)
def sum_bigwig_cmd(inbw):
    bw_sum, scale_factor = bwsum(inbw)
    print("bw_sum: {} | scale_factor: {}".format(bw_sum, scale_factor))


####################### Scale bigwig ##########################################
@cli.command("scale-bw", context_settings=CONTEXT_SETTINGS, help="Scale bigwig")
@click.option("--inbw", type=str, help="Path to input bigwig", required=True)
@click.option("--chromsizes", type=str, help="Path to chrom.sizes", required=True)
@click.option("--scalefactor", type=float, help="Scale factor", required=True)
@click.option("--outbw", type=str, help="Path to output bigwig", required=True)
def scale_bigwig_cmd(inbw, chromsizes, scalefactor, outbw):
    scale_bigwig(inbw, chromsizes, outbw, scalefactor)


####################### Normalize bigwig from hdf ##########################################
@cli.command(
    "normalize-bw-hdf",
    context_settings=CONTEXT_SETTINGS,
    help="Scale fragment length specific bigwig",
)
@click.option("--inbw", type=str, help="Path to input bigwig", required=True)
@click.option("--hdf", type=str, help="Path to HDF", required=True)
@click.option("--readlength", type=int, help="Fragment length to use", required=True)
@click.option("--outbw", type=str, help="Path to output bigwig", required=True)
def normalize_bw_hdf_cmd(inbw, hdf, readlength, outbw):
    normalize_bw_hdf(inbw, hdf, readlength, outbw)


######################### Plot periodicty dataframe ##############################
@cli.command(
    "plot-periodicity-df",
    context_settings=CONTEXT_SETTINGS,
    help="Plot periodicity heatmap for fragment length specific dataframe",
)
@click.option("--df", type=str, help="Path to tsv", required=True)
@click.option("--saveto", type=str, help="Output pdf/png", required=True)
def plot_periodicity_df_cmd(df, saveto):
    plot_periodicity_df(df, saveto)


######################### Shift bigwig ###########################################
@cli.command(
    "shift-bw",
    context_settings=CONTEXT_SETTINGS,
    help="Shift a bigwig by certain bases",
)
@click.option("--bw", type=str, help="Path to bigwig", required=True)
@click.option("--shift", type=int, help="Shift by these many bases", required=True)
@click.option("--outbw", type=str, help="Path to output bigwig", required=True)
def bwshift_cmd(bw, shift, outbw):
    bwshift(bw, shift, outbw)


####################### Create Loom files ######################################
@cli.command("loomify", context_settings=CONTEXT_SETTINGS, help="Create loom files")
@click.option("--root-dir", type=str, help="Root directory", required=True)
@click.option("--ribocop-dir", type=str, help="Pattern to search for", required=True)
@click.option(
    "--annotation", type=str, help="Path to RiboCop annotation file to", required=True
)
@click.option(
    "--out-dir", type=str, help="Path to RiboCop annotation file to", required=True
)
@click.option("--batch-size", type=int, help="Shift by these many bases", default=50)
def create_loom_files(root_dir, ribocop_dir, annotation, out_dir, batch_size):
    orf_tsv_list = glob.glob(
        "{}/*/{}/*_translating_ORFs.tsv".format(root_dir, ribocop_dir)
    )
    srx = [
        path_leaf(filepath).replace("_translating_ORFs.tsv", "")
        for filepath in orf_tsv_list
    ]
    srp = [filepath.split(ribocop_dir)[0].split("/")[-2] for filepath in orf_tsv_list]
    sample_list = list(zip(srp, srx, orf_tsv_list))
    write_loom_batches(sample_list, annotation, out_dir, batch_size=batch_size)


###################### orfs-seq function #########################################
@cli.command(
    "orfs-seq",
    context_settings=CONTEXT_SETTINGS,
    help="Generate sequence for ORFs in ribotricer's index",
)
@click.option(
    "--ribotricer_index",
    help=(
        "Path to the index file of ribotricer\n"
        "This file should be generated using ribotricer prepare-orfs"
    ),
    required=True,
)
@click.option("--fasta", help="Path to FASTA file", required=True)
@click.option("--saveto", help="Path to output file", required=True)
@click.option(
    "--offset_5p",
    help="Number of upstream bases to count(5')",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "--offset_3p",
    help="Number of downstream bases to count(3')",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "--protein", help="Output protein sequence instead of nucleotide", is_flag=True
)
def orf_seq_cmd(ribotricer_index, fasta, saveto, offset_5p, offset_3p, protein):
    if not os.path.isfile(ribotricer_index):
        sys.exit("Error: ribotricer index file not found")

    if not os.path.isfile(fasta):
        sys.exit("Error: fasta file not found")

    orf_seq(
        ribotricer_index,
        fasta,
        saveto,
        upstream_5p_offset=offset_5p,
        downstream_3p_offset=offset_3p,
        translate=protein,
    )
