import pandas as pd

from protrend.io import read_from_stack
from protrend.model import Effector
from protrend.transform.mix_ins import EffectorMixIn
from protrend.transform.regulondb.base import RegulonDBTransformer, regulondb_reader
from protrend.transform.transformations import drop_empty_string, drop_duplicates, create_input_value
from protrend.utils import SetList
from protrend.utils.processors import (rstrip, lstrip, apply_processors, remove_html_tags,
                                       parse_effector_name_regulondb)


class EffectorTransformer(EffectorMixIn, RegulonDBTransformer,
                          source='regulondb',
                          version='0.0.0',
                          node=Effector,
                          order=100,
                          register=True):
    default_transform_stack = {'effector': 'effector.txt'}
    columns = SetList(['protrend_id', 'name', 'kegg_compounds',
                       'effector_id', 'effector_name', 'category', 'effector_type', 'effector_note',
                       'effector_internal_comment', 'key_id_org'])
    read_columns = SetList(['effector_id', 'effector_name', 'category', 'effector_type', 'effector_note',
                            'effector_internal_comment', 'key_id_org'])

    @staticmethod
    def transform_effector(effector: pd.DataFrame):
        effector = effector.assign(name=effector['effector_name'].copy())

        effector = effector.dropna(subset=['name'])
        effector = drop_empty_string(effector, 'name')
        effector = drop_duplicates(df=effector, subset=['name'])

        effector = apply_processors(effector,
                                    effector_id=[rstrip, lstrip],
                                    name=[rstrip, lstrip, remove_html_tags, parse_effector_name_regulondb])

        effector = create_input_value(effector, 'name')
        return effector

    def transform(self):
        reader = regulondb_reader(skiprows=34, names=self.read_columns)
        effector = read_from_stack(stack=self.transform_stack, key='effector', columns=self.read_columns,
                                   reader=reader)

        effector = self.transform_effector(effector)
        annotated_effectors = self.annotate_effectors(effector)

        df = pd.merge(annotated_effectors, effector, on='input_value', suffixes=('_annotation', '_regulondb'))

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
