from typing import List

import pandas as pd

from protrend.io import read_txt, read_from_stack
from protrend.model.model import Effector
from protrend.transform.annotation import annotate_effectors
from protrend.transform.dto import EffectorDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, remove_html_tags, \
    parse_effector_name_regulondb
from protrend.transform.regulondb.base import RegulondbTransformer


class EffectorTransformer(RegulondbTransformer):
    default_node = Effector
    default_node_factors = ('name', )
    default_transform_stack = {'effector': 'effector.txt'}
    default_order = 100
    columns = {'protrend_id',
               'name', 'synonyms', 'mechanism', 'kegg_compounds',
               'effector_id', 'effector_name', 'category', 'effector_type', 'effector_note',
               'effector_internal_comment', 'key_id_org'}
    read_columns = {'effector_id', 'effector_name', 'category', 'effector_type', 'effector_note',
                    'effector_internal_comment', 'key_id_org'}

    def _transform_effector(self, effector: pd.DataFrame):
        effector = self.drop_duplicates(df=effector, subset=['name'], perfect_match=True, preserve_nan=True)
        effector = effector.dropna(subset=['name'])

        effector = apply_processors(effector, effector_id=[rstrip, lstrip],
                                    name=[rstrip, lstrip, remove_html_tags, parse_effector_name_regulondb])

        effector = self.create_input_value(effector, 'name')

        return effector

    @staticmethod
    def _transform_effectors(names: List[str]):
        dtos = [EffectorDTO(input_value=name) for name in names]
        annotate_effectors(dtos=dtos, names=names)

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        effector = read_from_stack(stack=self.transform_stack, file='effector', default_columns=self.read_columns,
                                   reader=read_txt, skiprows=34, names=self.read_columns)
        effector = self._transform_effector(effector)

        names = effector['input_value'].tolist()
        effectors = self._transform_effectors(names)

        df = pd.merge(effectors, effector, on='input_value')

        df = df.drop(columns=['input_value'])

        self._stack_transformed_nodes(df)

        return df
