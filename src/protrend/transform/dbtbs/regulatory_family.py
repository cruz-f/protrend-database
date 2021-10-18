import pandas as pd

from protrend.io import read_from_stack, read_json_lines
from protrend.model.model import RegulatoryFamily
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.transform.processors import (apply_processors, rstrip, lstrip, take_first)
from protrend.utils import SetList


class RegulatoryFamilyTransformer(DBTBSTransformer):
    default_node = RegulatoryFamily
    default_transform_stack = {'tf': 'TranscriptionFactor.json'}
    default_order = 100
    columns = SetList(['mechanism', 'name', 'description', 'rfam', 'protrend_id', 'tf', 'family'])
    read_columns = SetList(['name', 'family', 'domain', 'domain_description', 'description', 'url',
                            'type', 'comment', 'operon', 'subti_list', 'consensus_sequence'])

    def _transform_tf(self, tf: pd.DataFrame) -> pd.DataFrame:
        tf = tf.dropna(subset=['family'])
        tf = apply_processors(tf, family=[take_first, rstrip, lstrip])
        tf = self.drop_duplicates(df=tf, subset=['family'], perfect_match=True, preserve_nan=True)

        # drop not assigned
        tf_not_assigned_mask = tf['family'] == 'Not assigned'
        tf = tf[tf_not_assigned_mask]

        # drop other family
        tf_other_mask = tf['family'] == 'Other family'
        tf = tf[tf_other_mask]

        tf = tf.rename(columns={'name': 'tf'})
        tf = self.select_columns(tf, 'tf', 'family')

        tf['name'] = tf['family']

        # sigma factors
        tf_sigma_mask = tf['family'] == 'Sigma factors'
        tf.loc[tf_sigma_mask, 'mechanism'] = 'sigma factor'

        # tf mechanism
        tf_not_sigma_mask = tf['family'] != 'Sigma factors'
        tf.loc[tf_not_sigma_mask, 'mechanism'] = 'transcription factor'

        tf['rfam'] = None
        tf['description'] = None

        return tf

    def transform(self):
        tf = read_from_stack(stack=self.transform_stack, file='tf',
                             default_columns=self.read_columns, reader=read_json_lines)
        tf = self._transform_tf(tf)

        self._stack_transformed_nodes(tf)
        return tf
