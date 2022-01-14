import pandas as pd

from protrend.io import read, read_json_lines
from protrend.model import RegulatoryFamily, Regulator
from protrend.transform.dbtbs.base import DBTBSTransformer, DBTBSConnector
from protrend.transform.transformations import drop_empty_string, select_columns, group_by
from protrend.utils import SetList
from protrend.utils.processors import (apply_processors, rstrip, lstrip, take_first, to_list_nan, to_set_list)


class RegulatoryFamilyTransformer(DBTBSTransformer,
                                  source='dbtbs',
                                  version='0.0.4',
                                  node=RegulatoryFamily,
                                  order=100,
                                  register=True):
    columns = SetList(['protrend_id', 'name', 'mechanism', 'description',
                       'tf', 'family'])

    @staticmethod
    def transform_tf(tf: pd.DataFrame) -> pd.DataFrame:
        tf = apply_processors(tf, family=[to_list_nan, take_first, rstrip, lstrip])

        tf = tf.dropna(subset=['family'])
        tf = drop_empty_string(tf, 'family')

        # drop not assigned
        mask = tf['family'] != 'Not assigned'
        tf = tf[mask]

        # drop other family
        mask = tf['family'] != 'Other family'
        tf = tf[mask]

        tf = select_columns(tf, 'name', 'family')
        tf = tf.rename(columns={'name': 'tf', 'family': 'name'})

        tf = group_by(df=tf, column='name', aggregation={'tf': to_set_list})

        tf = tf.assign(mechanism=None, rfam=None, description=None)

        # sigma factors
        tf_sigma_mask = tf['name'] == 'Sigma factors'
        tf.loc[tf_sigma_mask, 'mechanism'] = 'sigma factor'

        # tf mechanism
        tf.loc[~tf_sigma_mask, 'mechanism'] = 'transcription factor'
        return tf

    def transform(self):
        tf = read(source=self.source, version=self.version, file='TranscriptionFactor.json', reader=read_json_lines,
                  default=pd.DataFrame(columns=['name', 'family', 'domain', 'domain_description', 'description', 'url',
                                                'type', 'consensus_sequence', 'comment', 'subti_list', 'gene', 'tfbs']))
        tf = self.transform_tf(tf)

        self.stack_transformed_nodes(tf)
        return tf


class RegulatoryFamilyToRegulatorConnector(DBTBSConnector,
                                           source='dbtbs',
                                           version='0.0.4',
                                           from_node=RegulatoryFamily,
                                           to_node=Regulator,
                                           register=True):

    def connect(self):
        source_df, target_df = self.transform_stacks(source='rfam',
                                                     target='regulator',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_on='tf',
                                                     target_on='dbtbs_name',
                                                     source_processors={'tf': [to_list_nan]},
                                                     target_processors={})
        source_df = source_df.explode('tf')
        source_df = apply_processors(source_df, tf=[rstrip, lstrip])

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='tf', target_on='dbtbs_name')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_connections(df)
