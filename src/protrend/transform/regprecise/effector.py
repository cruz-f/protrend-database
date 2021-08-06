from typing import List

import pandas as pd

from protrend.io.csv import read_csv
from protrend.io.json import read_json_lines
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

    def _transform_effector(self):

        file_path = self._transform_stack.get('effector')

        if not file_path:
            return pd.DataFrame(columns=['name', 'input_value'])

        df = read_json_lines(file_path)

        df = self.drop_duplicates(df=df, subset=['name'], perfect_match=True, preserve_nan=False)

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

    def transform(self):

        effector = self._transform_effector()

        names = list(effector['input_value'])

        effectors = self._transform_effectors(names)

        df = pd.merge(effectors, effector, on='input_value', suffixes=('_annotation', '_regprecise'))

        df['name'] = df['name_annotation'].astype(str) + df['name_regprecise'].astype(str)

        df = df.drop(['input_value', 'name_annotation', 'name_regprecise'], axis=1)

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df

    def _connect_to_source(self) -> pd.DataFrame:

        from_path = self._connect_stack.get('from')
        to_path = self._connect_stack.get('to')

        if not from_path:
            return pd.DataFrame()

        if not to_path:
            return pd.DataFrame()

        from_df = read_csv(from_path)
        from_identifiers = from_df['protrend_id'].tolist()

        size = len(from_identifiers)

        to_df = read_csv(to_path)
        to_df = to_df.query('name == regprecise')
        regprecise_id = to_df['protrend_id'].iloc[0]
        to_identifiers = [regprecise_id] * size

        kwargs = dict(url=from_df['url'].tolist(),
                      external_identifier=from_df['effector_id'].tolist(),
                      key=['effector_id'] * size)

        return self.make_connection(size=size,
                                    from_node=self.node,
                                    to_node=Source,
                                    from_identifiers=from_identifiers,
                                    to_identifiers=to_identifiers,
                                    kwargs=kwargs)

    def connect(self):

        source_connection = self._connect_to_source()
        df_name = f'connected_{self.node.node_name()}_{Source.node_name()}'
        self.stack_csv(df_name, source_connection)
