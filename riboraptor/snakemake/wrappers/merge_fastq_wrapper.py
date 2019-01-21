import os
from snakemake.shell import shell

input = (' ').join(snakemake.input.dynamic_input)

shell(r'''cat {input} > {snakemake.output}''')
