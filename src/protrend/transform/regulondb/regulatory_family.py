import pandas as pd

from protrend.io import read_from_stack, read_txt
from protrend.model.model import RegulatoryFamily
from protrend.transform.processors import (apply_processors, rstrip, lstrip)
from protrend.transform.regulondb.base import RegulondbTransformer
from protrend.utils import SetList


class RegulatoryFamilyTransformer(RegulondbTransformer):
    default_node = RegulatoryFamily
    default_transform_stack = {'tf': 'transcription_factor.txt'}
    default_order = 100
    columns = SetList(['transcription_factor_id', 'transcription_factor_name', 'site_length', 'symmetry',
                       'transcription_factor_family', 'tf_internal_comment', 'key_id_org',
                       'transcription_factor_note', 'connectivity_class', 'sensing_class', 'consensus_sequence',
                       'mechanism', 'name', 'protrend_id'])
    read_columns = SetList(['transcription_factor_id', 'transcription_factor_name', 'site_length', 'symmetry',
                            'transcription_factor_family', 'tf_internal_comment', 'key_id_org',
                            'transcription_factor_note', 'connectivity_class', 'sensing_class', 'consensus_sequence'])

    def _transform_tf(self, tf: pd.DataFrame) -> pd.DataFrame:
        tf = apply_processors(tf, transcription_factor_family=[rstrip, lstrip])
        tf = tf.dropna(subset=['transcription_factor_family'])
        tf = self.drop_duplicates(df=tf, subset=['transcription_factor_family'], perfect_match=True, preserve_nan=True)

        tf['name'] = tf['transcription_factor_family']
        tf['mechanism'] = 'transcription factor'

        return tf

    def transform(self):
        tf = read_from_stack(stack=self.transform_stack, file='tf',
                             default_columns=self.read_columns, reader=read_txt,
                             skiprows=38, names=self.read_columns)
        tf = self._transform_tf(tf)

        self._stack_transformed_nodes(tf)
        return tf
