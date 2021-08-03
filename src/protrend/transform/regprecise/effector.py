from typing import List

import pandas as pd

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
        return self._read_json_lines()

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

        names = tuple(effector['input_value'])

        effectors = self._transform_effectors(names)

        df = pd.merge(effectors, effector, on='input_value', suffixes=('_annotation', '_regprecise'))

        df['name'] = df['name_annotation'].astype(str) + df['name_regprecise'].astype(str)

        df = df.drop(['input_value', 'name_annotation', 'name_regprecise'], axis=1)

        n_rows, _ = df.shape
        df['source_db'] = ['regprecise'] * n_rows
        df['api_key'] = ['effector_id'] * n_rows

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df
