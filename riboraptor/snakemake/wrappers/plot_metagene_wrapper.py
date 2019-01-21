import os
from snakemake.shell import shell
if snakemake.wildcards.orientation == '5prime':
    RANGE = '-60:100'
else:
    RANGE = '-100:60'
if os.stat(str(snakemake.input)).st_size:
    shell(r'''riboraptor plot-metagene \
            --counts {snakemake.input} \
            --saveto {snakemake.output} \
            --positions {RANGE}''')
else:
    # Just touch the file
    shell(r'''touch {snakemake.output}''')
