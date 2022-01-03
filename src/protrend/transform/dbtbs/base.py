from abc import abstractmethod

import pandas as pd

from protrend.transform import Transformer, Connector


class DBTBSTransformer(Transformer, source='dbtbs', version='0.0.4', register=False):

    @abstractmethod
    def transform(self):
        pass

    def merge_annotations(self, annotations: pd.DataFrame, dbtbs: pd.DataFrame):
        df = pd.merge(annotations, dbtbs, on='input_value', suffixes=('_annotation', '_dbtbs'))

        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_dbtbs')

        df = df.assign(name_dbtbs_merge=df['name_dbtbs'].copy())
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_dbtbs_merge')

        df = self.merge_columns(df=df, column='genbank_accession',
                                left='genbank_accession_annotation', right='genbank_accession_dbtbs')

        df = self.merge_columns(df=df, column='uniprot_accession',
                                left='uniprot_accession_annotation', right='uniprot_accession_dbtbs')
        return df


class DBTBSConnector(Connector, source='dbtbs', version='0.0.4', register=False):

    @abstractmethod
    def connect(self):
        pass
