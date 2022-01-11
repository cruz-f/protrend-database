from abc import abstractmethod

import pandas as pd

from protrend.transform import Transformer, Connector
from protrend.transform.transformations import merge_columns, select_columns


class RegPreciseTransformer(Transformer, source='regprecise', version='0.0.0', register=False):

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

    @abstractmethod
    def transform(self):
        pass


class RegPreciseConnector(Connector, source='regprecise', version='0.0.0', register=False):

    @abstractmethod
    def connect(self):
        pass
