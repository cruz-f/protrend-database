from typing import List

import pandas as pd

from protrend.model.model import Source
from protrend.transform.annotation.effector import annotate_effectors
from protrend.transform.dto import EffectorDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, nan_to_str
from protrend.transform.regprecise.settings import EffectorSettings
from protrend.transform.transformer import Transformer


class EffectorTransformer(Transformer):

    def __init__(self, settings: EffectorSettings = None):

        if not settings:
            settings = EffectorSettings()

        super().__init__(settings)

    def read(self, **kwargs):
        files = self._read_json_lines()

        return super(EffectorTransformer, self).read(**files)

    @staticmethod
    def _transform_effector(df):
        df = df.drop_duplicates(subset=['name'])

        apply_processors(rstrip, lstrip, df=df, col='name')

        df['input_value'] = df['name']

        return df

    @staticmethod
    def _transform_effectors(names: List[str]):

        dtos = [EffectorDTO(input_value=name) for name in names]
        annotate_effectors(dtos=dtos, names=names)

        effectors = pd.DataFrame([dto.to_dict() for dto in dtos])

        if effectors.empty:
            effectors = pd.DataFrame(columns=['input_value', 'name'])

        apply_processors(nan_to_str, df=effectors, col='name')

        return effectors

    def transform(self, **kwargs):

        effector = kwargs.get('effector', pd.DataFrame(columns=['name']))
        effector = self._transform_effector(effector)

        names = list(effector['input_value'])

        effectors = self._transform_effectors(names)

        df = pd.merge(effectors, effector, on='input_value', suffixes=('_annotation', '_regprecise'))

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
                      external_identifier=list(df['effector_id']),
                      name=['effector_id'] * size)

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
