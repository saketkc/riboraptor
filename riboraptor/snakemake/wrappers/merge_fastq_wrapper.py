from snakemake.shell import shell

if len(snakemake.input.dynamic_input) > 1:
    shell(r'''
          cat {snakemake.input.dynamic_input} > {snakemake.output}
          ''')
else:
    shell(r'''
          ln -s {snakemake.input.dynamic_input} {snakemake.output}
          ''')
