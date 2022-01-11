import pandas as pd

from protrend.io import read_from_stack
from protrend.model import RegulatoryFamily, Regulator
from protrend.transform.regulondb.base import RegulonDBTransformer, RegulonDBConnector, regulondb_reader
from protrend.transform.transformations import drop_empty_string, drop_duplicates
from protrend.utils import SetList
from protrend.utils.processors import (apply_processors, rstrip, lstrip)


class RegulatoryFamilyTransformer(RegulonDBTransformer,
                                  source='regulondb',
                                  version='0.0.0',
                                  node=RegulatoryFamily,
                                  order=100,
                                  register=True):
    default_transform_stack = {'tf': 'transcription_factor.txt'}
    columns = SetList(['protrend_id', 'name', 'mechanism', 'description',
                       'transcription_factor_id', 'transcription_factor_name', 'site_length', 'symmetry',
                       'transcription_factor_family', 'tf_internal_comment', 'key_id_org',
                       'transcription_factor_note', 'connectivity_class', 'sensing_class', 'consensus_sequence'])
    read_columns = SetList(['transcription_factor_id', 'transcription_factor_name', 'site_length', 'symmetry',
                            'transcription_factor_family', 'tf_internal_comment', 'key_id_org',
                            'transcription_factor_note', 'connectivity_class', 'sensing_class', 'consensus_sequence'])

    @staticmethod
    def transform_tf(tf: pd.DataFrame) -> pd.DataFrame:
        tf = tf.assign(name=tf['transcription_factor_family'].copy(), mechanism='transcription factor',
                       description=None)

        tf = apply_processors(tf, name=[rstrip, lstrip])
        tf = tf.dropna(subset=['name'])
        tf = drop_empty_string(tf, 'name')
        tf = drop_duplicates(df=tf, subset=['name'])

        # removing families with _like and _
        name = tf['name'].str.replace('_like', '')
        name = name.str.replace('_', ' ')
        tf = tf.assign(name=name)

        return tf

    def transform(self):
        tf_reader = regulondb_reader(skiprows=38, names=self.read_columns)
        tf = read_from_stack(stack=self.transform_stack, key='tf',
                             columns=self.read_columns, reader=tf_reader)
        tf = self.transform_tf(tf)

        self.stack_transformed_nodes(tf)
        return tf


class RegulatoryFamilyToRegulatorConnector(RegulonDBConnector,
                                           source='regulondb',
                                           version='0.0.0',
                                           from_node=RegulatoryFamily,
                                           to_node=Regulator,
                                           register=True):
    default_connect_stack = {'rfam': 'integrated_regulatoryfamily.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        df = self.create_connection(source='rfam', target='regulator',
                                    source_on='transcription_factor_id', target_on='transcription_factor_id')
        self.stack_json(df)
