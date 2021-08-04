from typing import List

import pandas as pd

from protrend.model.model import Source
from protrend.transform.annotation.organism import annotate_organisms
from protrend.transform.dto import OrganismDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, nan_to_str
from protrend.transform.regprecise.settings import OrganismSettings
from protrend.transform.transformer import Transformer


class OrganismTransformer(Transformer):

    def __init__(self, settings: OrganismSettings = None):

        if not settings:
            settings = OrganismSettings()

        super().__init__(settings)

    def read(self, **kwargs):
        files = self._read_json_lines()

        return super(OrganismTransformer, self).read(**files)

    @staticmethod
    def _transform_genome(df):

        df = df.drop_duplicates(subset=['name'])

        apply_processors(rstrip, lstrip, df=df, col='name')

        df['input_value'] = df['name']

        return df

    @staticmethod
    def _transform_organisms(names: List[str]):

        dtos = [OrganismDTO(input_value=name) for name in names]
        annotate_organisms(dtos=dtos, names=names)

        organisms = pd.DataFrame([dto.to_dict() for dto in dtos])

        if organisms.empty:
            organisms = pd.DataFrame(columns=['input_value', 'name'])

        apply_processors(nan_to_str, df=organisms, col='name')

        return organisms

    def transform(self, **kwargs):

        genome = kwargs.get('genome', pd.DataFrame(columns=['name']))
        genome = self._transform_genome(genome)

        names = list(genome['input_value'])

        organisms = self._transform_organisms(names)

        df = pd.merge(organisms, genome, on='input_value', suffixes=('_annotation', '_regprecise'))

        df['name'] = df['name_annotation'].astype(str) + df['name_regprecise'].astype(str)

        df = df.drop(['input_value', 'name_annotation', 'name_regprecise'], axis=1)

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df

    def _connect_to_source(self, df: pd.DataFrame) -> pd.DataFrame:

        size, _ = df.shape

        from_identifiers = list(df[self.node.identifying_property])
        to_identifiers = ['regprecise'] * size

        kwargs = dict(url=list(df['url']),
                      external_identifier=list(df['genome_id']),
                      name=['genome_id'] * size)

        return self.make_connection(size=size,
                                    from_node=self.node,
                                    to_node=Source,
                                    from_property=self.node.identifying_property,
                                    to_property='name',
                                    from_identifiers=from_identifiers,
                                    to_identifiers=to_identifiers,
                                    **kwargs)

    def connect(self, df: pd.DataFrame):

        source_connection = self._connect_to_source(df)
        df_name = f'connected_{self.node.node_name()}_{Source.node_name()}'
        self.stack_csv(df_name, source_connection)