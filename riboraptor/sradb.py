"""Helper functions for parsing SRAmetadb.sqlite file"""

import re
import sqlite3
import pandas as pd
"""
Tables

1. Study
> dbListFields(sra_con,"study")
"study_ID"    "study_alias"    "study_accession"    "study_title"    "study_type"
"study_abstract"    "broker_name"    "center_name"    "center_project_name"
"study_description"    "related_studies"    "primary_study"    "sra_link"
"study_url_link"    "xref_link"    "study_entrez_link"    "ddbj_link"    "ena_link"


2. Sample
> dbListFields(sra_con,"sample")
"sample_ID"    "sample_alias"    "sample_accession"    "broker_name"
"center_name"    "taxon_id"    "scientific_name"    "common_name"
"anonymized_name"    "individual_name"    "description"    "sra_link"
"sample_url_link"    "xref_link"    "sample_entrez_link"   "ddbj_link"
"ena_link"    "sample_attribute"    "submission_accession" "sradb_updated"

3. Experiment
> dbListFields(sra_con,"experiment")
"experiment_ID"    "bamFile"    "fastqFTP"    "experiment_alias"
"experiment_accession"    "broker_name"    "center_name"    "title"
"study_name"    "study_accession"    "design_description"    "sample_name"
"sample_accession"    "sample_member"    "library_name"    "library_strategy"
"library_source"    "library_selection"    "library_layout"    "targeted_loci"
"library_construction_protocol" "spot_length"    "adapter_spec"    "read_spec"
"platform"    "instrument_model"    "platform_parameters"    "sequence_space"
"base_caller"    "quality_scorer"    "number_of_levels"    "multiplier"
"qtype"    "sra_link"    "experiment_url_link"    "xref_link"
"experiment_entrez_link"    "ddbj_link"    "ena_link"    "experiment_attribute"
"submission_accession"    "sradb_updated"

4. Run
> dbListFields(sra_con,"run")
"run_ID"    "bamFile"    "run_alias"    "run_accession"
"broker_name"    "instrument_name"    "run_date"    "run_file"
"run_center"    "total_data_blocks"    "experiment_accession"    "experiment_name"
"sra_link"    "run_url_link"    "xref_link"    "run_entrez_link"     "ddbj_link"
"ena_link"    "run_attribute"    "submission_accession" "sradb_updated"


5. Submission

> dbListFields(sra_con,"submission")
"submission_ID"    "submission_alias"    "submission_accession"   "submission_comment"
"files"    "broker_name"    "center_name"    "lab_name"
"submission_date"    "sra_link"    "submission_url_link"    "xref_link"
"submission_entrez_link"   "ddbj_link"    "ena_link"    "submission_attribute"
"sradb_updated"

"""


def _extract_first_field(data):
    """Extract first field from a list of fields"""
    return list(next(zip(*data)))


class SRAdb(object):
    def __init__(self, sqlite_file):
        self.sqlite_file = sqlite_file
        self.open()
        self.cursor = self.db.cursor()
        self.valid_in_acc_type = [
            'SRA', 'ERA', 'DRA', 'SRP', 'ERP', 'DRP', 'SRS', 'ERS', 'DRS',
            'SRX', 'ERX', 'DRX', 'SRR', 'ERR', 'DRR'
        ]
        self.valid_in_type = {
            'SRA': 'submission',
            'ERA': 'submission',
            'DRA': 'submission',
            'SRP': 'study',
            'ERP': 'study',
            'DRP': 'study',
            'SRS': 'sample',
            'ERS': 'sample',
            'DRS': 'sample',
            'SRX': 'experiment',
            'ERX': 'experiment',
            'DRX': 'experiment',
            'SRR': 'run',
            'ERR': 'run',
            'DRR': 'run'
        }

    def open(self):
        self.db = sqlite3.connect(self.sqlite_file)
        self.db.text_factory = str

    def close(self):
        self.db.close()

    def list_tables(self):
        """List all tables in the sqlite"""
        results = self.cursor.execute(
            'SELECT name FROM sqlite_master WHERE type="table";').fetchall()
        return _extract_first_field(results)

    def list_fields(self, table):
        "List all fields in a given table"
        results = self.cursor.execute('SELECT * FROM {}'.format(table))
        return _extract_first_field(results.description)

    def desc_table(self, table):
        results = self.cursor.execute(
            'PRAGMA table_info("{}")'.format(table)).fetchall()
        columns = ['cid', 'name', 'dtype', 'notnull', 'dflt_value', 'pk']
        data = []
        for result in results:
            data.append(list(map(lambda x: str(x), result)))
        df = pd.DataFrame(data, columns=columns)
        return df

    def get_query(self, query):
        results = self.cursor.execute(query).fetchall()
        column_names = list(map(lambda x: x[0], self.cursor.description))
        results = [dict(zip(column_names, result)) for result in results]
        return pd.DataFrame(results)

    def get_row_count(self, table):
        """Get row counts for a table"""
        return self.cursor.execute(
            'SELECT max(rowid) FROM {}'.format(table)).fetchone()[0]

    def get_table_counts(self):
        tables = self.list_tables()
        results = dict(
            [(table, self.get_row_count(table)) for table in tables])
        return pd.DataFrame.from_dict(
            results, orient='index', columns=['count'])

    def sra_convert(self,
                    acc,
                    out_type=[
                        'study_accession',
                        'experiment_accession',
                        'run_accession',
                        'taxon_id',
                        'library_selection',
                        'library_layout',
                        'library_strategy',
                        'library_source',
                        'library_name',
                        'bases',
                        'spots',
                        'adapter_spec',
                    ]):
        in_acc_type = re.sub('\\d+$', '', acc).upper()
        if in_acc_type not in self.valid_in_acc_type:
            raise ValueError('{} not a valid input type'.format(in_acc_type))

        in_type = self.valid_in_type[in_acc_type]
        out_type = [x for x in out_type if x != in_type]

        select_type = [in_type + '_accession'] + out_type
        select_type_sql = (',').join(select_type)
        sql = "SELECT DISTINCT " + select_type_sql + " FROM sra_ft WHERE sra_ft MATCH '" + acc + "';"
        df = self.get_query(sql)
        df['avg_read_length'] = df['bases'] / df['spots']
        df['spots'] = df['spots'].astype(int)
        df['bases'] = df['bases'].astype(int)
        df['taxon_id'] = df['taxon_id'].fillna(0).astype(int)
        df = df.sort_values(by=[
            'taxon_id', 'avg_read_length', 'run_accession',
            'experiment_accession', 'library_selection'
        ])
        df = df[out_type + ['avg_read_length']].reset_index(drop=True)
        return df

    def search_experiment(self, srx):
        """Search for a SRX/GSM id in the experiments"""
        if 'GSM' in srx:
            results = self.cursor.execute(
                'select * from EXPERIMENT where experiment_alias = "{}"'.
                format(srx)).fetchall()
        else:
            results = self.cursor.execute(
                'select * from EXPERIMENT where experiment_accession = "{}"'.
                format(srx)).fetchall()
        assert len(results) == 1, 'Got multiple hits'
        results = results[0]
        column_names = list(map(lambda x: x[0], self.cursor.description))
        results = dict(zip(column_names, results))
        return results
