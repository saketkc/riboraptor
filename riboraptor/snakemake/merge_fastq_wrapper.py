from snakemake.shell import shell

if len(snakemake.input) > 1:
    shell(r'''
          cat {snakemake.input} > {snakemake.output}
          ''')
else:
    shell(r'''
          ln -s {snakemake.input} {snakemake.output}
          ''')
