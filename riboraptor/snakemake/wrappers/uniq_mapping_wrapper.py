import tempfile
from snakemake.shell import shell
with tempfile.TemporaryDirectory(dir=snakemake.params.tmp_dir) as temp_dir:
    shell(r'''samtools view -b -q 255 \
          {snakemake.input} -o {snakemake.output}.temp \
          && samtools sort -@ {snakemake.threads} \
          {snakemake.output}.temp -o {snakemake.output} \
          -T {temp_dir}/{snakemake.wildcards.sample}_sort \
          && rm -rf {snakemake.output}.temp \
          && samtools index {snakemake.output}
          ''')
