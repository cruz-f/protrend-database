from abc import abstractmethod

import pandas as pd

from protrend.transform import Transformer, Connector
from protrend.transform.transformations import merge_columns, locus_tag_curation


class DBTBSTransformer(Transformer, register=False):

    @abstractmethod
    def transform(self):
        pass

    @staticmethod
    def merge_annotations(annotations: pd.DataFrame, dbtbs: pd.DataFrame):
        df = pd.merge(annotations, dbtbs, on='input_value', suffixes=('_annotation', '_dbtbs'))

        df = merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_dbtbs')

        df = merge_columns(df=df, column='name', left='name_annotation', right='name_dbtbs')

        df = merge_columns(df=df, column='genbank_accession',
                           left='genbank_accession_annotation', right='genbank_accession_dbtbs')

        df = merge_columns(df=df, column='uniprot_accession',
                           left='uniprot_accession_annotation', right='uniprot_accession_dbtbs')

        df = locus_tag_curation(df)
        return df


class DBTBSConnector(Connector, register=False):

    @abstractmethod
    def connect(self):
        pass
