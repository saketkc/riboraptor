import os
from textwrap import dedent
import base64


def generate_base64_images(list_of_images):
    text = '<tr>'
    for m in list_of_images:
        m = str(m)
        if not os.path.isfile(m) or not os.stat(str(m)).st_size:
            continue
        with open(m, 'rb') as f:
            encoded = base64.b64encode(f.read())
        img_tag = '<img alt="" src="data:image/png;base64,{0}" height="400px">'.format(
            str((encoded), 'utf-8'))
        text += dedent('''
                       <td>
                       <figure>
                       {img_tag}
                       </figure>
                       </td>''').format(img_tag=img_tag)
    text += '</tr>'

    return text


text = dedent('''<html>
              <center>
              <h1> Summary </h1>
              <table style="width:100%">
              <tr>
                <th>Fragment Length</th>
                <th>Periodicity v/s Fragment Length</th>
                <th>Metagene (All fragments)</th>
              </tr>
              ''')

for metagene, fragment_length, periodicity in zip(
        snakemake.input['metagene'], snakemake.input['fragment_length'],
        snakemake.input['periodicity']):
    text += generate_base64_images([fragment_length, periodicity, metagene])

text += dedent('''</table></center>
               </html>''')

with open(snakemake.output.html, 'w') as f:
    f.write(text)
