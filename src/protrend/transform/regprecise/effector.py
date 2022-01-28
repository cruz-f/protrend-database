import pandas as pd

from protrend.io import read_json_lines, read
from protrend.model import Effector
from protrend.transform.mix_ins import EffectorMixIn
from protrend.transform.regprecise.base import RegPreciseTransformer
from protrend.transform.transformations import drop_empty_string, drop_duplicates, create_input_value
from protrend.utils import SetList
from protrend.utils.processors import rstrip, lstrip, apply_processors, to_int_str


class EffectorTransformer(EffectorMixIn, RegPreciseTransformer,
                          source='regprecise',
                          version='0.0.0',
                          node=Effector,
                          order=100,
                          register=True):
    columns = SetList(['protrend_id', 'name', 'kegg_compounds',
                       'effector_id', 'url', 'regulog'])

    @staticmethod
    def transform_effector(effector: pd.DataFrame):
        # noinspection DuplicatedCode
        effector = effector.dropna(subset=['name'])
        effector = drop_empty_string(effector, 'name')
        effector = drop_duplicates(df=effector, subset=['name'])

        effector = apply_processors(effector, effector_id=to_int_str, name=[rstrip, lstrip])

        effector = create_input_value(effector, 'name')
        return effector

    def transform(self):
        effector = read(source=self.source, version=self.version,
                        file='Effector.json', reader=read_json_lines,
                        default=pd.DataFrame(columns=['effector_id', 'name', 'url', 'regulog']))

        effectors = self.transform_effector(effector)
        annotated_effectors = self.annotate_effectors(effectors)

        df = self.merge_annotations_by_name(annotated_effectors, effectors)

        self.stack_transformed_nodes(df)
        return df
