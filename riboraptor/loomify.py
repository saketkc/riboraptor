"""Module to store RiboCop output in loom format"""

from ast import literal_eval
import os
from tqdm import tqdm

from .helpers import mkdir_p

import numpy as np
import pandas as pd
import loompy


def _create_index_for_annotation_row(row):
    """Parse the row from annotation file to get chrom and positions.

    Parameters
    ----------
    row: pd.DataFrame row
         Row of annotation dataframe

    Returns
    -------
    index: array_like
           Array of 1-based positions on chromosome
    """
    coordinates = row["coordinate"].split(",")
    index = []
    for coordinate in coordinates:
        coordinate_start, coordinate_end = [int(x) for x in coordinate.split("-")]
        index += list(range(coordinate_start, coordinate_end + 1))
    return index


def batch(iterable, n=50):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx : min(ndx + n, l)]


def read_batch_tsv(sample_list):
    dfs = [
        pd.read_table(
            tsv,
            usecols=["ORF_ID", "profile"],
            dtype={"ORF_ID": "str", "profile": "str"},
            memory_map=True,
        )
        for srp, srx, tsv in sample_list
    ]
    srps = [srp for srp, srx, tsv in sample_list]
    srxs = [srx for srp, srx, tsv in sample_list]
    col_attrs = {"study": srps, "experiment": srxs}
    return dfs, col_attrs


def write_loom_file(loom_file_path, matrix, col_attrs, row_attrs=None):
    """Create/Update loom file.

    Parameters
    ----------
    filepath: string
              Path to write loomfile
    matrix: array_like
            Data matrix
    col_attrs: dict
               A dict of lists with same length as the columns in matrix
    row_attrs: dict
               A dict of lists with same length as the rows in matrix

    """
    # Check if loomfile exists
    if os.path.isfile(loom_file_path):
        # Read its columns attributes
        # to see if the new column_attrs are same as previous
        with loompy.connect(loom_file_path) as ds:
            assert sorted(list(col_attrs.keys())) == sorted(list(ds.col_attrs.keys()))
            # Check if nothing is being added again
            for key in col_attrs.keys():
                common_columns = list(
                    set(col_attrs[key]).intersection(set(ds.col_attrs[key]))
                )
                if len(common_columns) > 0:
                    # Delete common columns
                    for common_column in common_columns:
                        # Find this index of key
                        index = col_attrs[key].index(common_column)
                        # Delete this column from matrix
                        matrix = np.delete(matrix, index, 1)
                        # Remove it from the list
                        col_attrs[key].remove(common_column)
                        assert matrix.shape[1] == len(col_attrs[key])

        # Add columns
        with loompy.connect(loom_file_path) as dsout:
            for key in col_attrs.keys():
                col_attrs[key] = np.array(col_attrs[key])
                assert matrix.shape[1] == len(col_attrs[key])
            dsout.add_columns(matrix, col_attrs=col_attrs)
    else:
        for key in col_attrs.keys():
            col_attrs[key] = np.array(col_attrs[key])
        loompy.create(loom_file_path, matrix, row_attrs=row_attrs, col_attrs=col_attrs)


def get_row_attrs_from_orf(annotation, orf_id):
    """Get row attributes from annoation file for given orf_id.

    Parameters
    ----------
    annotation: pd.DataFrame
                Annotation dataframe

    orf_id: str
            ORF_ID
    """
    orf_row = annotation.loc[annotation.ORF_ID == orf_id].iloc[0]
    orf_index = _create_index_for_annotation_row(orf_row)
    chrom = orf_row["chrom"]
    row_attrs = ["{}_{}".format(chrom, pos) for pos in orf_index]
    row_attrs = {"chr_pos": np.array(row_attrs)}
    return row_attrs


def write_loom_for_dfs(loom_file_path, list_of_dfs, col_attrs, row_attrs, orf_id):
    """Write loom file for a list of dataframes.

    Parameters
    ----------
    loom_file_path: str
                    Path to output loomfile
    list_of_dfs:  list
                  List of dataframes
    col_attrs: dict
               A dict of lists with same length as the columns in matrix
    row_attrs: dict
               A dict of lists with same length as the rows in matrix
    orf_id: str
            ORF_ID
    """
    dfs_subset = [df[df.ORF_ID == orf_id].iloc[0] for df in list_of_dfs]
    profile_stacked = [literal_eval(df.profile) for df in dfs_subset]
    matrix = np.array(profile_stacked).T
    write_loom_file(loom_file_path, matrix, col_attrs, row_attrs)


def write_loom_batches(sample_list, annotation_filepath, out_root_dir, batch_size=50):
    print("Reading annotation ... ")
    annotation = pd.read_table(annotation_filepath)
    print("Done!")

    ORF_IDS = annotation.ORF_ID.tolist()
    TX_IDS = annotation.transcript_id.tolist()
    with tqdm(total=len(sample_list) // batch_size) as pbar:
        for batch_sample in batch(sample_list, batch_size):
            # Read all the tsvs in this batch
            print("Reading batch tsv ... ")
            list_of_dfs, col_attrs = read_batch_tsv(batch_sample)
            print("Done! ")
            # For each ORF in all tsvs
            # write a loom file
            # organized by `transcript_id/ORF_ID.loom`
            for tx_id, orf_id in zip(TX_IDS, ORF_IDS):
                row_attrs = get_row_attrs_from_orf(annotation, orf_id)
                out_dir = os.path.join(out_root_dir, tx_id)
                mkdir_p(out_dir)
                loom_file_path = os.path.join(out_dir, "{}.loom".format(orf_id))
                write_loom_for_dfs(
                    loom_file_path, list_of_dfs, col_attrs, row_attrs, orf_id
                )
            del list_of_dfs
            del col_attrs
            pbar.update()
