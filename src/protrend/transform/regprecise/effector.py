from typing import List

import pandas as pd

from protrend.io.utils import read_from_stack
from protrend.transform.annotation.effector import annotate_effectors
from protrend.transform.connector import DefaultConnector
from protrend.transform.dto import EffectorDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors
from protrend.transform.regprecise.settings import EffectorSettings, EffectorToSource
from protrend.transform.regprecise.source import SourceTransformer
from protrend.transform.transformer import DefaultTransformer


class EffectorTransformer(DefaultTransformer):
    default_settings = EffectorSettings
    columns = {'protrend_id', 'effector_id', 'name', 'url', 'regulog', 'mechanism', 'synonyms', 'kegg_compounds'}
    read_columns = {'effector_id', 'name', 'url', 'regulog'}

    def _transform_effector(self, effector: pd.DataFrame):
        df = self.drop_duplicates(df=effector, subset=['name'], perfect_match=True, preserve_nan=False)

        apply_processors(rstrip, lstrip, df=df, col='name')

        df['input_value'] = df['name']

        return df

    @staticmethod
    def _transform_effectors(names: List[str]):
        dtos = [EffectorDTO(input_value=name) for name in names]
        annotate_effectors(dtos=dtos, names=names)

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        effector = read_from_stack(tl=self, file='effector', json=True, default_columns=self.read_columns)
        effector = self._transform_effector(effector)

        names = list(effector['input_value'])

        effectors = self._transform_effectors(names)

        df = pd.merge(effectors, effector, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise', fill='')

        df = df.drop(['input_value'], axis=1)

        if df.empty:
            df = self.make_empty_frame()

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df


class EffectorToSourceConnector(DefaultConnector):
    default_settings = EffectorToSource

    def connect(self):
        effector = read_from_stack(tl=self, file='effector', json=False, default_columns=EffectorTransformer.columns)
        source = read_from_stack(tl=self, file='source', json=False, default_columns=SourceTransformer.columns)

        from_identifiers = effector['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=effector['url'].tolist(),
                      external_identifier=effector['effector_id'].tolist(),
                      key=['effector_id'] * size)

        df = self.make_connection(size=size,
                                  from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_csv(df)
