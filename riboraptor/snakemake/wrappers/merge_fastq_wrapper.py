import os
from snakemake.shell import shell

if len(snakemake.input.dynamic_input) > 1:
    shell(r'''
          cat {snakemake.input.dynamic_input} > {snakemake.output}
          ''')
else:
    source = os.path.abspath(snakemake.input.dynamic_input)
    destination = os.path.dirname(os.path.abspath(snakemake.output))
    shell(r'''
          ln -s {source} {destination}
          ''')
