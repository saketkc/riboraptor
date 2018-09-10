import os
from textwrap import dedent
from riboraptor.helpers import path_leaf
import base64
sample = snakemake.wildcards.sample


def generate_base64_images(list_of_images):
    text = ''
    for m in list_of_images:
        m = str(m)
        if not os.path.isfile(m) or not os.stat(str(m)).st_size:
            continue
        fig_name = path_leaf(m).replace('.png', '')
        fragment_length = path_leaf(os.path.dirname(m))
        with open(m, 'rb') as f:
            encoded = base64.b64encode(f.read())
        img_tag = '<img alt="" src="data:image/png;base64,{0}" height="400px">'.format(
            str((encoded), 'utf-8'))
        text += dedent('''
        <figure>
        {img_tag}
        <figcaption> {fig_name} </figcaption>
        </figure>''').format(
            img_tag=img_tag,
            fig_name=fig_name + ' | Read Length: {}'.format(fragment_length))
    return text


text = dedent('''<html>
              <center>
              <h1> Sample: {sample} </h1>
              <h1> Fragment Length Distribution </h1>
              '''.format(sample=sample))
text += generate_base64_images([snakemake.input.fragment_length])

text += dedent('''<h1> Metagene Plot (combined)</h1>
               ''')
text += generate_base64_images([snakemake.input.metagene])

text += dedent('''<h1> Metagene Plots (5' positive strand)</h1>
               ''')
text += generate_base64_images(snakemake.input.prime5_pos)
text += dedent('''<h1> Metagene Plots (3' positive strand)</h1>
               ''')
text += generate_base64_images(snakemake.input.prime3_pos)

text += dedent('''<h1> Metagene Plots (5' negative strand)</h1>
               ''')
text += generate_base64_images(snakemake.input.prime5_neg)

text += dedent('''<h1> Metagene Plots (3' negative strand)</h1>
               ''')
text += generate_base64_images(snakemake.input.prime3_neg)

text += dedent('''<h1> Metagene Plots (5' both strands combined)</h1>
               ''')
text += generate_base64_images(snakemake.input.prime5_combined)

text += dedent('''<h1> Metagene Plots (3' both strands combined)</h1>
               ''')
text += generate_base64_images(snakemake.input.prime3_combined)

text += dedent('''</center>
               </html>''')

with open(snakemake.output.html, 'w') as f:
    f.write(text)
