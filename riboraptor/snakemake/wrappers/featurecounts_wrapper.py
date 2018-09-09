import os
import h5py
from snakemake.shell import shell
protocols = []
for f in snakemake.input.hdfs:
    hdf = h5py.File(f, 'r')
    protocol = hdf.attrs['protocol']
    protocols.append(protocol)
    hdf.close()
assert len(set(protocols)) == 1
protocol = protocols[0]
if protocol == 'forward':
    count_strat = '-s 1'
elif protocol == 'unstranded':
    count_strat = '-s 2'
else:
    count_strat = ''
bams = sorted(snakemake.input.bams)
shell(
    r'''featureCounts {count_strat} -a {snakemake.params.annotation} -o {snakemake.output} -t exon -g gene_id -Q 4 -T {snakemake.threads} {bams}'''
)
