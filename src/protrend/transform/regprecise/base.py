from abc import abstractmethod

import pandas as pd

from protrend.io import read, read_json_lines
from protrend.transform import Transformer, Connector
from protrend.transform.transformations import merge_columns, select_columns


class RegPreciseTransformer(Transformer, register=False):

    @staticmethod
    def merge_annotations_by_name(annotation: pd.DataFrame, original: pd.DataFrame) -> pd.DataFrame:
        df = pd.merge(annotation, original, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        df = df.drop(columns=['input_value'])
        return df

    @staticmethod
    def transform_organism(organism: pd.DataFrame) -> pd.DataFrame:
        # see regulator transform organism method
        organism = select_columns(organism, 'protrend_id', 'name', 'ncbi_taxonomy', 'genome_id')

        taxonomy = organism['ncbi_taxonomy'].fillna(organism['name'])
        organism = organism.assign(ncbi_taxonomy=taxonomy)
        organism = organism.drop(columns=['name'])
        organism = organism.rename(columns={'protrend_id': 'organism', 'genome_id': 'genome'})
        return organism

    def transform_rfmas(self) -> pd.DataFrame:
        tf_family_columns = ['tffamily_id', 'name', 'url', 'description', 'pubmed', 'regulog']
        tf_columns = ['collection_id', 'name', 'url', 'description', 'pubmed', 'regulog']
        rna_columns = ['riboswitch_id', 'name', 'url', 'description', 'pubmed', 'rfam', 'regulog']

        tf_family = read(source=self.source, version=self.version,
                         file='TranscriptionFactorFamily.json', reader=read_json_lines,
                         default=pd.DataFrame(columns=tf_family_columns))

        tf = read(source=self.source, version=self.version,
                  file='TranscriptionFactor.json', reader=read_json_lines,
                  default=pd.DataFrame(columns=tf_columns))

        rna = read(source=self.source, version=self.version,
                   file='RNAFamily.json', reader=read_json_lines,
                   default=pd.DataFrame(columns=rna_columns))

        tf_family = self.transform_tf_family(tf_family)
        tf = self.transform_tf(tf)
        rna = self.transform_rna(rna)

        df = pd.concat([tf_family, tf, rna])
        df = df.reset_index(drop=True)
        return df

    def transform_tf_family(self, tf_family: pd.DataFrame) -> pd.DataFrame:
        pass

    def transform_tf(self, tf: pd.DataFrame) -> pd.DataFrame:
        pass

    def transform_rna(self, rna: pd.DataFrame) -> pd.DataFrame:
        pass

    @abstractmethod
    def transform(self):
        pass


class RegPreciseConnector(Connector, register=False):

    @abstractmethod
    def connect(self):
        pass
