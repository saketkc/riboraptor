import os
import tempfile
from snakemake.shell import shell

if len(snakemake.input) > 1:
    with tempfile.TemporaryDirectory(dir=snakemake.params.tmp_dir) as temp_dir:
        cmd = ' -in '.join(snakemake.input)
        shell(r'''bamtools merge -in {cmd} -out {snakemake.output}.unsorted \
              && samtools sort -@ {snakemake.threads} \
              -T {temp_dir}/{snakemake.wildcards.sample}_merge_bam \
              -o {snakemake.output} {snakemake.output}.unsorted \
              && samtools index {snakemake.output} \
              && yes | rm -rf {snakemake.output}.unsorted''')
elif len(snakemake.input) == 1:
    source = os.path.abspath(str(snakemake.input[0]))
    destination = os.path.abspath(str(snakemake.output))
    shell('''cp {source} {destination} && cp {source}.bai {destination}.bai''')
