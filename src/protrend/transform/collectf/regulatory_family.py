import pandas as pd

from protrend.io import read_from_stack, read_json_lines
from protrend.model.model import RegulatoryFamily
from protrend.transform.collectf.base import CollectfTransformer
from protrend.transform.processors import apply_processors, remove_white_space, rstrip, lstrip, \
    remove_multiple_white_space, parse_collectf_description


class RegulatoryFamilyTransformer(CollectfTransformer):
    default_node = RegulatoryFamily
    default_node_factors = ('name',)
    default_transform_stack = {'tf': 'TranscriptionFactor.json'}
    default_order = 100
    columns = {'protrend_id',
               'mechanism'
               'name', 'family', 'description', 'regulon'}
    read_columns = {'name', 'family', 'description', 'regulon'}

    def _transform_tf(self, tf: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=tf, subset=['name'], perfect_match=True, preserve_nan=True)

        df = apply_processors(df, name=remove_white_space,
                              description=[parse_collectf_description, remove_multiple_white_space, rstrip, lstrip])

        df['mechanism'] = 'transcription factor'

        return df

    def transform(self):
        tf = read_from_stack(stack=self.transform_stack, file='tf',
                             default_columns=self.read_columns, reader=read_json_lines)
        tf = self._transform_tf(tf)

        self._stack_transformed_nodes(tf)
        return tf
