import pandas as pd

from protrend.io import read_from_stack, read_txt, read_json_frame
from protrend.model import RegulatoryFamily, Regulator
from protrend.utils.processors import (apply_processors, rstrip, lstrip)
from protrend.transform.regulondb.base import RegulondbTransformer, RegulondbConnector
from protrend.transform.regulondb.regulator import RegulatorTransformer
from protrend.utils import SetList


class RegulatoryFamilyTransformer(RegulondbTransformer,
                                  source='regulondb',
                                  version='0.0.0',
                                  node=RegulatoryFamily,
                                  order=100,
                                  register=True):
    default_transform_stack = {'tf': 'transcription_factor.txt'}
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
        tf = self.drop_duplicates(df=tf, subset=['transcription_factor_family'], perfect_match=True)

        tf['name'] = tf['transcription_factor_family']
        tf['mechanism'] = 'transcription factor'

        return tf

    def transform(self):
        tf = read_from_stack(stack=self.transform_stack, file='tf',
                             default_columns=self.read_columns, reader=read_txt,
                             skiprows=38, names=self.read_columns)
        tf = self._transform_tf(tf)

        self.stack_transformed_nodes(tf)
        return tf


class RegulatorToRegulatoryFamilyConnector(RegulondbConnector,
                                           source='regulondb',
                                           version='0.0.0',
                                           from_node=Regulator,
                                           to_node=RegulatoryFamily,
                                           register=True):
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'rfam': 'integrated_regulatoryfamily.json'}

    def connect(self):
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = regulator[['protrend_id', 'transcription_factor_id']]
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        rfam = read_from_stack(stack=self._connect_stack, file='rfam',
                               default_columns=RegulatoryFamilyTransformer.columns, reader=read_json_frame)
        rfam = rfam[['protrend_id', 'transcription_factor_id']]
        rfam = rfam.rename(columns={'protrend_id': 'rfam_protrend_id'})

        df = pd.merge(regulator, rfam, on='transcription_factor_id')
        df = df.drop_duplicates(subset=['regulator_protrend_id', 'rfam_protrend_id'])

        from_identifiers = df['regulator_protrend_id'].tolist()
        to_identifiers = df['rfam_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
