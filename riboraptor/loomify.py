"""Module to store RiboCop output in loom format"""

from ast import literal_eval
import os
from joblib import Parallel, delayed
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
    assert len(srps) > 0
    assert len(srxs) > 0
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
            for key in col_attrs.keys():
                assert len(col_attrs[key]) > 0
            # Check if nothing is being added again
            items_to_delete = []
            index_to_delete = []
            for key in col_attrs.keys():
                common_columns = list(
                    set(col_attrs[key]).intersection(set(ds.col_attrs[key]))
                )
                if len(common_columns) > 0:
                    items_to_delete = []
                    # Delete common columns
                    for common_column in common_columns:
                        # Find this index of key
                        index = col_attrs[key].index(common_column)
                        index_to_delete.append(index)
                        items_to_delete.append(common_column)

                    print("len(common_column): {}".format(len(common_columns)))
                    print("common_column: {}".format(common_columns))

                    for item in items_to_delete:
                        # Remove it from the list
                        col_attrs[key].remove(item)
                        print("common_Col_lene: {}".format(len(col_attrs[key])))
        print("matrix.shape: {}".format(matrix.shape))
        print("len(index_to_delete): {}".format(len(index_to_delete)))
        print("len(items_to_delete): {}".format(len(items_to_delete)))
        for index in set(index_to_delete):
            # Delete this column from matrix
            matrix = np.delete(matrix, index, 1)
        print("matrix.shape: {}".format(matrix.shape))
        for key in col_attrs.keys():
            print("len(col_attrs[{}]): {}".format(key, len(col_attrs[key])))

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
    print("##########col_attrs:{}".format(col_attrs))
    dfs_subset = [df[df.ORF_ID == orf_id].iloc[0] for df in list_of_dfs]
    profile_stacked = [literal_eval(df.profile) for df in dfs_subset]
    profile_ncol = len(profile_stacked)
    profile_nrow = len(profile_stacked[0])
    matrix = np.array(profile_stacked).T
    # print('#######matrix########: {}'.format(matrix))
    matrix = matrix.reshape(profile_nrow, profile_ncol)
    # print('nrow: {}'.format(profile_nrow))
    # print('ncol: {}'.format(profile_ncol))
    # print('shape: {}'.format(matrix.shape))
    write_loom_file(loom_file_path, matrix, col_attrs, row_attrs)


def _run_parallel_writer(data):
    orf_id = data["orf_id"]
    row_attrs = data["orf_attrs"]
    loom_file_path = data["loom_file_path"]
    list_of_dfs = data["list_of_dfs"]
    col_attrs = data["col_attrs"]
    mkdir_p(os.path.dirname(loom_file_path))
    write_loom_for_dfs(loom_file_path, list_of_dfs, col_attrs.copy(), row_attrs, orf_id)


def write_loom_batches(sample_list, annotation_filepath, out_root_dir, batch_size=50):
    print("Reading annotation ... ")
    annotation = pd.read_table(annotation_filepath)
    print("Done!")

    ORF_IDS = annotation.ORF_ID.tolist()
    TX_IDS = annotation.transcript_id.tolist()
    ROW_ATTRS = [get_row_attrs_from_orf(annotation, orf_id) for orf_id in ORF_IDS]
    OUT_DIRS = [os.path.join(out_root_dir, tx_id) for tx_id in TX_IDS]
    LOOM_FILE_PATHS = [
        os.path.join(out_dir, "{}.loom".format(orf_id))
        for out_dir, orf_id in zip(OUT_DIRS, ORF_IDS)
    ]
    DATA = []
    print("Processing dict data")
    for tx_id, orf_id, row_attrs, out_dir, loom_file_path in zip(
        tqdm(TX_IDS), ORF_IDS, ROW_ATTRS, OUT_DIRS, LOOM_FILE_PATHS
    ):
        data_dict = {
            "orf_id": orf_id,
            "row_attrs": row_attrs,
            "loom_file_path": loom_file_path,
        }
        DATA.append(data_dict)
    del TX_IDS
    del LOOM_FILE_PATHS
    del OUT_DIRS
    del ROW_ATTRS
    print("Done!")
    with tqdm(total=len(sample_list) // batch_size) as pbar:
        for batch_sample in batch(sample_list, batch_size):
            # Read all the tsvs in this batch
            print("Reading batch tsv ... ")
            list_of_dfs, col_attrs = read_batch_tsv(batch_sample)
            print("############ col_attrs:{}".format(col_attrs))
            print("Done! ")
            # For each ORF in all tsvs
            # write a loom file
            # organized by `transcript_id/ORF_ID.loom`
            print("Changing dict")
            for data in tqdm(DATA):
                data["list_of_dfs"] = list_of_dfs
                data["col_attrs"] = col_attrs
            print("Done")

            Parallel(n_jobs=16)(delayed(_run_parallel_writer)(data) for data in DATA)
            """
            for tx_id, orf_id, row_attr, out_dir, loom_file_path in zip(TX_IDS, ORF_IDS, ROW_ATTRS, OUT_DIRS, LOOM_FILE_PATHS):
                mkdir_p(out_dir)
                write_loom_for_dfs(
                    loom_file_path, list_of_dfs, col_attrs.copy(), row_attrs, orf_id
                )
            """
            del list_of_dfs
            del col_attrs
            pbar.update()
