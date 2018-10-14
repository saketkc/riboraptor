import os
import pandas as pd
import numpy as np
import pybedtools
import tempfile
from riboraptor.helpers import path_leaf
from snakemake.shell import shell


def total_genome_size(chrom_sizes_file):
    """Return sum total of chromosome sizes"""
    df = pd.read_table(chrom_sizes_file, names=['chrom', 'sizes'])
    total = df['sizes'].sum()
    return total


def get_align_intro_params(intron_bed_file):
    df = pybedtools.BedTool(intron_bed_file).to_dataframe()
    lengths = df['end'] - df['start']

    ## Based on small genomes. See https://groups.google.com/forum/#!topic/rna-star/hQeHTBbkc0c
    alignintronNmin = max(4, lengths.min())
    alignintronNmax = lengths.max()
    return alignintronNmin, alignintronNmax


OUT_PREFIX = os.path.splitext(snakemake.output[0])[0]
TMP_DIR_SAMPLE = path_leaf(OUT_PREFIX)
STAR_LOGS_DIR = os.path.dirname(snakemake.output.starlogs)

ALIGN_INTRON_Nmin, ALIGN_INTRON_Nmax = get_align_intro_params(
    snakemake.params.intron_bed)
TOTAL_GENOME_SIZE = total_genome_size(snakemake.params.chrom_sizes)
SA_INDEX_Nbases = int(np.floor(min(14, np.log2(TOTAL_GENOME_SIZE) / 2.0 - 1)))

with tempfile.TemporaryDirectory(dir=snakemake.params.tmp_dir) as temp_dir:
    shell(r'''
            STAR --runThreadN {snakemake.threads}\
                --genomeDir {snakemake.input.index}\
                --outFilterMismatchNmax 2\
                --alignIntronMin {ALIGN_INTRON_Nmin}\
                --alignIntronMax {ALIGN_INTRON_Nmax}\
                --outFileNamePrefix {OUT_PREFIX}\
                --readFilesIn {snakemake.input.R1}\
                --readFilesCommand zcat\
                --quantMode TranscriptomeSAM GeneCounts\
                --outSAMtype BAM Unsorted\
                --outTmpDir {temp_dir}/{TMP_DIR_SAMPLE}_tmp\
                --outFilterType BySJout\
                --outFilterMatchNmin 16\
                --seedSearchStartLmax 15\
                --winAnchorMultimapNmax 200\
                && samtools sort -@ {snakemake.threads} {OUT_PREFIX}Aligned.out.bam -o {snakemake.output.bam} -T {temp_dir}/{TMP_DIR_SAMPLE}_sort\
                && mv {OUT_PREFIX}Aligned.toTranscriptome.out.bam {snakemake.output.txbam}\
                && samtools index {snakemake.output.bam}\
                && mv {OUT_PREFIX}ReadsPerGene.out.tab {snakemake.output.counts}\
                && mv {OUT_PREFIX}Log.final.out {OUT_PREFIX}Log.out {OUT_PREFIX}SJ.out.tab\
                {OUT_PREFIX}Log.progress.out {STAR_LOGS_DIR}
            ''')
