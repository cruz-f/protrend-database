import pandas as pd

from protrend.io import read_from_stack, read_json_lines, read_json_frame
from protrend.model.model import RegulatoryFamily, Regulator
from protrend.transform.dbtbs.base import DBTBSTransformer, DBTBSConnector
from protrend.transform.dbtbs.regulator import RegulatorTransformer
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
        tf = tf[~tf_not_assigned_mask]

        # drop other family
        tf_other_mask = tf['family'] == 'Other family'
        tf = tf[~tf_other_mask]

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


class RegulatorToRegulatoryFamilyConnector(DBTBSConnector):
    default_from_node = Regulator
    default_to_node = RegulatoryFamily
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'rfam': 'integrated_regulatoryfamily.json'}

    def connect(self):
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = regulator[['protrend_id', 'name_dbtbs']]
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        rfam = read_from_stack(stack=self._connect_stack, file='rfam',
                               default_columns=RegulatoryFamilyTransformer.columns, reader=read_json_frame)
        rfam = rfam[['protrend_id', 'tf']]
        rfam = rfam.rename(columns={'protrend_id': 'rfam_protrend_id'})

        df = pd.merge(regulator, rfam, right_on='name_dbtbs', left_on='tf')
        df = df.drop_duplicates(subset=['regulator_protrend_id', 'rfam_protrend_id'])

        from_identifiers = df['regulator_protrend_id'].tolist()
        to_identifiers = df['rfam_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)