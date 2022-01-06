import pandas as pd

from protrend.io import read_json_lines, read_from_stack
from protrend.model import Effector
from protrend.transform.regprecise.base import RegPreciseTransformer
from protrend.utils import SetList
from protrend.utils.processors import rstrip, lstrip, apply_processors, to_int_str


class EffectorTransformer(RegPreciseTransformer,
                          source='regprecise',
                          version='0.0.0',
                          node=Effector,
                          order=100,
                          register=True):
    default_transform_stack = {'effector': 'Effector.json'}
    columns = SetList(['protrend_id', 'name', 'kegg_compounds',
                       'effector_id', 'url', 'regulog'])
    read_columns = SetList(['effector_id', 'name', 'url', 'regulog'])

    def transform_effector(self, effector: pd.DataFrame):
        # noinspection DuplicatedCode
        effector = effector.dropna(subset=['name'])
        effector = self.drop_empty_string(effector, 'name')
        effector = self.drop_duplicates(df=effector, subset=['name'])

        effector = apply_processors(effector, effector_id=to_int_str, name=[rstrip, lstrip])

        effector = self.create_input_value(effector, 'name')
        return effector

    def transform(self):
        effector = read_from_stack(stack=self.transform_stack, key='effector', columns=self.read_columns,
                                   reader=read_json_lines)

        effectors = self.transform_effector(effector)
        annotated_effectors = self.annotate_effectors(effectors)

        df = self.merge_annotations_by_name(annotated_effectors, effectors)

        self.stack_transformed_nodes(df)
        return df
