import tempfile
from snakemake.shell import shell
with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
    shell(r'''samtools view -b -q 255 {input} -o {output}.temp \
          && samtools sort -@ {threads} {output}.temp -o {output} \
          -T {temp_dir}/{wildcards.sample}_sort \
          && rm -rf {output}.temp \
          && samtools index {output}
          ''')
