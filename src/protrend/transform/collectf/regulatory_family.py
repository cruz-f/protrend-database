import pandas as pd

from protrend.io import read_from_stack, read_json_lines, read_json_frame
from protrend.model import RegulatoryFamily, Regulator
from protrend.transform.collectf.base import CollectfTransformer, CollectfConnector
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.utils.processors import (apply_processors, remove_white_space, rstrip, lstrip,
                                       remove_multiple_white_space, parse_collectf_description, to_list)
from protrend.utils import SetList


class RegulatoryFamilyTransformer(CollectfTransformer,
                                  source='collectf',
                                  version='0.0.1',
                                  node=RegulatoryFamily,
                                  order=100,
                                  register=True):
    default_transform_stack = {'tf': 'TranscriptionFactor.json'}
    columns = SetList(['name', 'family', 'description', 'regulon', 'mechanism', 'protrend_id'])
    read_columns = SetList(['name', 'family', 'description', 'regulon'])

    def _transform_tf(self, tf: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=tf, subset=['name'], perfect_match=True)

        df = apply_processors(df,
                              name=remove_white_space,
                              description=[parse_collectf_description, remove_multiple_white_space, rstrip, lstrip])

        df['mechanism'] = 'transcription factor'

        return df

    def transform(self):
        tf = read_from_stack(stack=self.transform_stack, file='tf',
                             default_columns=self.read_columns, reader=read_json_lines)
        tf = self._transform_tf(tf)

        self.stack_transformed_nodes(tf)
        return tf


class RegulatoryFamilyToRegulatorConnector(CollectfConnector,
                                           source='collectf',
                                           version='0.0.1',
                                           from_node=RegulatoryFamily,
                                           to_node=Regulator,
                                           register=True):
    default_connect_stack = {'rfam': 'integrated_regulatoryfamily.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        regulator = read_from_stack(stack=self.connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        rfam = read_from_stack(stack=self.connect_stack, file='rfam',
                               default_columns=RegulatoryFamilyTransformer.columns, reader=read_json_frame)
        rfam = rfam.rename(columns={'protrend_id': 'rfam_protrend_id'})
        rfam = apply_processors(rfam, regulon=to_list)
        rfam = rfam.explode(column='regulon')

        df = pd.merge(regulator, rfam, left_on='uniprot_accession', right_on='regulon')

        from_identifiers = df['rfam_protrend_id'].tolist()
        to_identifiers = df['regulator_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)
        self.stack_json(df)
