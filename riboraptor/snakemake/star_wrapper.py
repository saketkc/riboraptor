import os
import pandas as pd
import numpy as np
import tempfile
from riboraptor.helpers import path_leaf
from snakemake.shell import shell


def total_genome_size(chrom_sizes_file):
    """Return sum total of chromosome sizes"""
    df = pd.read_table(chrom_sizes_file, names=['chrom', 'sizes'])
    total = df['sizes'].sum()
    return total


def get_align_intro_params(intron_bed_file):
    df = pd.read_table(
        intron_bed_file,
        names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    lengths = df['end'] - df['start']

    ## Based on small genomes. See https://groups.google.com/forum/#!topic/rna-star/hQeHTBbkc0c
    alignintronNmin = max(4, lengths.min())
    alignintronNmax = lengths.max()
    return alignintronNmin, alignintronNmax


OUT_PREFIX = os.path.splitext(snakemake.output[0])[0]
TMP_DIR_SAMPLE = path_leaf(OUT_PREFIX)
START_LOGS_DIR = os.dirname(snakemake.output.starlogs)

ALIGN_INTRON_Nmin, ALIGN_INTRON_Nmax = get_align_intro_params(
    snakemake.params.intron_bed)
TOTAL_GENOME_SIZE = total_genome_size(snakemake.params.chrom_sizes)
SA_INDEX_Nbases = int(np.floor(min(14, np.log2(TOTAL_GENOME_SIZE) / 2.0 - 1)))

with tempfile.TemporaryDirectory(dir=snakemake.params.tmp_dir) as temp_dir:
    shell(r'''
            STAR --runThreadN {snakemake.threads}'\
                --genomeDir {snakemake.input.index}\
                --outFilterMismatchNmax 2\
                --alignIntronMin {ALIGN_INTRON_Nmin}\
                --alignIntronMax {ALIGN_INTRON_Nmax}\
                --outFileNamePrefix {out_prefix}\
                --readFilesIn {snakemake.input.R1}\
                --readFilesCommand zcat\
                --quantMode TranscriptomeSAM GeneCounts\
                --outSAMtype BAM Unsorted\
                --outTmpDir {temp_dir}/{TMP_DIR_SAMPLE}_tmp\
                --outFilterType BySJout\
                --outFilterMatchNmin 16\
                --seedSearchStartLmax 15\
                --winAnchorMultimapNmax 200\
                && samtools sort -@ {snakemake.threads} {snakemake.params.prefix}Aligned.out.bam -o {snakemake.output.bam} -T {temp_dir}/{TMP_DIR_SAMPLE}_sort\
                && mv {snakemake.params.prefix}Aligned.toTranscriptome.out.bam {snakemake.output.txbam}\
                && samtools index {snakemake.output.bam}\
                && mv {snakemake.params.prefix}ReadsPerGene.out.tab {snakemake.output.counts}\
                && mv {params.prefix}Log.final.out {params.prefix}Log.out {params.prefix}SJ.out.tab\
                {params.prefix}Log.progress.out {STAR_LOGS_DIR}
            ''')
