import os
from textwrap import dedent
from riboraptor.helpers import path_leaf
import base64


def generate_base64_images(list_of_images):
    text = ''
    for m in list_of_images:
        m = str(m)
        if not os.path.isfile(m) or not os.stat(str(m)).st_size:
            continue
        fig_name = path_leaf(m).replace('.png', '')
        with open(m, 'rb') as f:
            encoded = base64.b64encode(f.read())
        img_tag = '<img alt="" src="data:image/png;base64,{0}" height="400px">'.format(
            str((encoded), 'utf-8'))
        text += dedent('''
        <figure>
        {img_tag}
        <figcaption> {fig_name} </figcaption>
        </figure>''').format(
            img_tag=img_tag, fig_name=fig_name)
    return text


text = dedent('''<html>
              <center>
              <h1> QC plots </h1>
              ''')
text = generate_base64_images(snakemake.input.dynamic_input)
text += dedent('''
               </center>
               </html>
               ''')
with open(snakemake.output.html, 'w') as f:
    f.write(text)
