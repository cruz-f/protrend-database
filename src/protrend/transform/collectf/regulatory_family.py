import pandas as pd

from protrend.io import read_json_lines, read
from protrend.model import RegulatoryFamily, Regulator
from protrend.transform.collectf.base import CollecTFTransformer, CollecTFConnector
from protrend.transform.transformations import select_columns, group_by
from protrend.utils import SetList
from protrend.utils.constants import TRANSCRIPTION_FACTOR
from protrend.utils.processors import apply_processors, rstrip, lstrip, flatten_set_list_nan, to_list_nan


class RegulatoryFamilyTransformer(CollecTFTransformer,
                                  source='collectf',
                                  version='0.0.1',
                                  node=RegulatoryFamily,
                                  order=100,
                                  register=True):
    columns = SetList(['protrend_id', 'name', 'mechanism', 'description',
                       'regulon'])

    @staticmethod
    def transform_tf(tf: pd.DataFrame) -> pd.DataFrame:
        tf = apply_processors(tf, family=[rstrip, lstrip], regulon=to_list_nan)
        tf = select_columns(tf, 'family', 'regulon')

        tf = group_by(df=tf, column='family', aggregation={}, default=flatten_set_list_nan)
        tf = tf.rename(columns={'family': 'name'})
        tf = tf.assign(mechanism=TRANSCRIPTION_FACTOR, description=None)
        return tf

    def transform(self):
        tf = read(source=self.source, version=self.version, file='TranscriptionFactor.json', reader=read_json_lines,
                  default=pd.DataFrame(columns=['name', 'family', 'description', 'regulon']))
        tf = self.transform_tf(tf)

        self.stack_transformed_nodes(tf)
        return tf


class RegulatoryFamilyToRegulatorConnector(CollecTFConnector,
                                           source='collectf',
                                           version='0.0.1',
                                           from_node=RegulatoryFamily,
                                           to_node=Regulator,
                                           register=True):

    def connect(self):
        source_df, target_df = self.transform_stacks(source='rfam',
                                                     target='regulator',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_on='regulon',
                                                     target_on='uniprot_accession',
                                                     source_processors={'regulon': [to_list_nan]},
                                                     target_processors={})
        source_df = source_df.explode('regulon')
        source_df = apply_processors(source_df, regulon=[rstrip, lstrip])

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='regulon', target_on='uniprot_accession')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_connections(df)
