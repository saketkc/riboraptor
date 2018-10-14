import os
import h5py
from snakemake.shell import shell
from collections import defaultdict
import pandas as pd

protocols = defaultdict(list)
for hdf, bam in zip(snakemake.input['hdfs'], snakemake.input['bams']):
    hdf = h5py.File(hdf, 'r')
    protocol = hdf.attrs['protocol']
    protocols[protocol].append(bam)
    hdf.close()
outputs = []
for protocol, bams in protocols.items():
    if protocol == 'forward':
        count_strat = '-s 1'
    elif protocol == 'unstranded':
        count_strat = '-s 2'
    else:
        count_strat = ''
    bams = sorted(bams)
    output = os.path.abspath(str(snakemake.output)) + '-' + protocol
    outputs.append(output)
    shell(
        r'''featureCounts {count_strat} -a {snakemake.params.annotation} -o {output} -t exon -g gene_id -Q 4 -T {snakemake.threads} {bams}'''
    )
df = pd.read_table(outputs[0], skiprows=[0])
if len(outputs) > 1:
    for f in outputs[1:]:
        temp_df = pd.read_table(f, skiprows=[0])
        temp_df = temp_df.drop(
            columns=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'])
        df = pd.concat([df, temp_df], axis=1)

df.to_csv(str(snakemake.output), sep='\t', index=False, header=True)
