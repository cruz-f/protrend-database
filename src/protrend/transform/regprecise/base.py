from abc import abstractmethod

import pandas as pd

from protrend.transform import Transformer, Connector


class RegPreciseTransformer(Transformer, source='regprecise', version='0.0.0', register=False):

    def merge_annotations_by_name(self, annotation: pd.DataFrame, original: pd.DataFrame) -> pd.DataFrame:
        df = pd.merge(annotation, original, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        df = df.drop(columns=['input_value'])
        return df

    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        # see regulator transform organism method
        organism = self.select_columns(organism, 'protrend_id', 'name', 'ncbi_taxonomy', 'genome_id')

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
