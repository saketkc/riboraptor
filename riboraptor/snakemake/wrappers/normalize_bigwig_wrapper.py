import os
from snakemake.shell import shell

if os.stat(str(snakemake.input.bw)).st_size:
    shell(r'''riboraptor normalize-bw-hdf \
          --inbw {snakemake.input.bw} \
          --hdf {snakemake.input.hdf} \
          --readlength {snakemake.wildcards.fragment_length} \
          --outbw {snakemake.output}
          ''')
else:
    # Just touch the file
    shell(r'''touch {snakemake.output}''')
