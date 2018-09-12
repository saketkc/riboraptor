import os
from riboraptor.helpers import path_leaf
from riboraptor.helpers import mkdir_p
from riboraptor.coherence import naive_periodicity
import pandas as pd
from collections import defaultdict

periodicities = defaultdict(list)
for f in snakemake.input:
    fragment_length = path_leaf(os.path.dirname(f))
    sample_name = path_leaf(os.path.dirname(os.path.dirname(f)))
    if not os.path.isfile(str(f)) or not os.stat(str(f)).st_size:
        continue
    counts = pd.read_table(f)
    counts = pd.Series(
        counts['count'].tolist(), index=counts['position'].tolist())
    periodicity = naive_periodicity(counts)
    periodicities[int(fragment_length)].append(periodicity)

df = pd.DataFrame.from_dict(
    periodicities, orient='index', columns=[sample_name]).sort_index()
mkdir_p(os.path.dirname(str(snakemake.output)))
df.to_csv(str(snakemake.output), sep='\t', index=True, header=True)
