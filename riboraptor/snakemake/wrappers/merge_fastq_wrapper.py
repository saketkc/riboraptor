import os
from snakemake.shell import shell

if len(snakemake.input.dynamic_input) > 1:
    shell(r'''
          cat {snakemake.input.dynamic_input} > {snakemake.output}
          ''')
else:
    source = os.path.abspath(str(snakemake.input.dynamic_input[0]))
    destination = os.path.abspath(str(snakemake.output))
    shell(r'''cp {source} {destination}''')
