import pandas as pd

from protrend.io.utils import read_rfam, read_regulator
from protrend.model import RegulatoryFamily, Regulator
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.transformations import select_columns, drop_empty_string, drop_duplicates
from protrend.utils import SetList
from protrend.utils.constants import TRANSCRIPTION_FACTOR, SMALL_RNA
from protrend.utils.processors import (remove_regprecise_more, remove_multiple_white_space,
                                       rstrip, lstrip, remove_pubmed, apply_processors, to_int_str)


class RegulatoryFamilyTransformer(RegPreciseTransformer,
                                  source='regprecise',
                                  version='0.0.0',
                                  node=RegulatoryFamily,
                                  order=100,
                                  register=True):
    columns = SetList(['protrend_id', 'name', 'mechanism', 'rfam', 'description',
                       'tffamily_id', 'collection_id', 'riboswitch_id', 'url'])

    def transform_tf_family(self, tf_family: pd.DataFrame) -> pd.DataFrame:
        tf_family = select_columns(tf_family, 'tffamily_id', 'name', 'url', 'description', 'pubmed')

        # noinspection DuplicatedCode
        tf_family = tf_family.dropna(subset=['name'])
        tf_family = drop_empty_string(tf_family, 'name')
        tf_family = drop_duplicates(df=tf_family, subset=['name'])

        tf_family = apply_processors(tf_family,
                                     name=[rstrip, lstrip],
                                     description=[remove_regprecise_more, remove_pubmed, remove_multiple_white_space,
                                                  rstrip, lstrip])

        tf_family = tf_family.assign(mechanism=TRANSCRIPTION_FACTOR)
        return tf_family

    def transform_tf(self, tf: pd.DataFrame) -> pd.DataFrame:
        tf = select_columns(tf, 'collection_id', 'name', 'url', 'description', 'pubmed')

        # noinspection DuplicatedCode
        tf = tf.dropna(subset=['name'])
        tf = drop_empty_string(tf, 'name')
        tf = drop_duplicates(df=tf, subset=['name'])

        tf = apply_processors(tf,
                              name=[rstrip, lstrip],
                              description=[remove_regprecise_more, remove_pubmed, remove_multiple_white_space,
                                           rstrip, lstrip])

        tf = tf.assign(mechanism=TRANSCRIPTION_FACTOR)
        return tf

    def transform_rna(self, rna: pd.DataFrame) -> pd.DataFrame:
        rna = select_columns(rna, 'riboswitch_id', 'name', 'url', 'description', 'rfam', 'pubmed')

        # noinspection DuplicatedCode
        rna = rna.dropna(subset=['name'])
        rna = drop_empty_string(rna, 'name')
        rna = drop_duplicates(df=rna, subset=['name'])

        rna = apply_processors(rna,
                               name=[rstrip, lstrip],
                               rfam=[rstrip, lstrip],
                               description=[remove_regprecise_more, remove_pubmed, remove_multiple_white_space, rstrip,
                                            lstrip])

        rna = rna.assign(mechanism=SMALL_RNA)
        return rna

    def transform(self):
        df = self.transform_rfmas()

        # clean the other regulatory family
        mask = df['name'] != '[Other]'
        df = df[mask]

        df = drop_duplicates(df, subset=['name'])

        df = apply_processors(df, tffamily_id=to_int_str, riboswitch_id=to_int_str, collection_id=to_int_str)

        self.stack_transformed_nodes(df)
        return df


class RegulatoryFamilyToRegulatorConnector(RegPreciseConnector,
                                           source='regprecise',
                                           version='0.0.0',
                                           from_node=RegulatoryFamily,
                                           to_node=Regulator,
                                           register=True):

    def connect(self):
        source = read_rfam(source=self.source, version=self.version, columns=RegulatoryFamilyTransformer.columns)
        target = read_regulator(source=self.source, version=self.version, columns=RegulatorTransformer.columns)

        tfs_chain = [('tffamily_id', 'tf_family'),
                     ('collection_id', 'transcription_factor'),
                     ('riboswitch_id', 'rna_family')]

        from_identifiers = []
        to_identifiers = []
        for source_key, target_key in tfs_chain:
            source_df = source.dropna(subset=[source_key])

            target_df = target.explode(target_key)
            target_df = target_df.dropna(subset=[target_key])

            merged = pd.merge(source_df, target_df, left_on=source_key, right_on=target_key,
                              suffixes=('_rfam', '_regulator'))

            merged = merged.dropna(subset=['protrend_id_rfam', 'protrend_id_regulator'])
            merged = merged.drop_duplicates(subset=['protrend_id_rfam', 'protrend_id_regulator'])

            from_identifiers.extend(merged['protrend_id_rfam'].to_list())
            to_identifiers.extend(merged['protrend_id_regulator'].to_list())

        df = self.connection_frame(source_ids=from_identifiers, target_ids=to_identifiers)
        self.stack_connections(df)
