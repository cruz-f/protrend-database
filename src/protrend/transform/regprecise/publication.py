import pandas as pd

from protrend.model import Publication, RegulatoryFamily
from protrend.transform.mix_ins import PublicationMixIn
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.transformations import (drop_empty_string, select_columns, group_by, create_input_value,
                                                merge_columns)
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, to_int_str, to_list_nan, to_set_list, flatten_set_list_nan


class PublicationTransformer(PublicationMixIn, RegPreciseTransformer,
                             source='regprecise',
                             version='0.0.0',
                             node=Publication,
                             order=100,
                             register=True):
    columns = SetList(['protrend_id', 'pmid', 'doi', 'title', 'author', 'year',
                       'tffamily_id', 'collection_id', 'riboswitch_id', 'pubmed'])

    @staticmethod
    def _transform_rfams(rfam: pd.DataFrame):
        rfam = rfam.assign(pmid=rfam['pubmed'].copy())

        rfam = apply_processors(rfam, pmid=to_list_nan)
        rfam = rfam.explode('pmid')

        rfam = rfam.dropna(subset=['pmid'])
        rfam = drop_empty_string(rfam, 'pmid')

        return rfam

    def transform_tf_family(self, tf_family: pd.DataFrame) -> pd.DataFrame:
        tf_family = select_columns(tf_family, 'tffamily_id', 'pubmed')
        return self._transform_rfams(tf_family)

    def transform_tf(self, tf: pd.DataFrame) -> pd.DataFrame:
        tf = select_columns(tf, 'collection_id', 'pubmed')
        return self._transform_rfams(tf)

    def transform_rna(self, rna: pd.DataFrame) -> pd.DataFrame:
        rna = select_columns(rna, 'riboswitch_id', 'pubmed')
        return self._transform_rfams(rna)

    def transform(self):
        publications = self.transform_rfmas()

        aggregation = {'tffamily_id': to_set_list, 'collection_id': to_set_list, 'riboswitch_id': to_set_list,
                       'pubmed': flatten_set_list_nan}
        publications = group_by(publications, column='pmid', aggregation=aggregation)
        publications = create_input_value(publications, col='pmid')

        annotated_publications = self.annotate_publications(publications)

        df = pd.merge(annotated_publications, publications, on='input_value', suffixes=('_annotation', '_regprecise'))

        # merge pmid
        df = merge_columns(df=df, column='pmid', left='pmid_annotation',
                           right='pmid_regprecise')

        df = apply_processors(df, pmid=to_int_str, year=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df


class PublicationToRegulatoryFamilyConnector(RegPreciseConnector,
                                             source='regprecise',
                                             version='0.0.0',
                                             from_node=Publication,
                                             to_node=RegulatoryFamily,
                                             register=True):

    def connect(self):
        source_df, target_df = self.transform_stacks(source='publication',
                                                     target='rfam',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_on='pubmed',
                                                     target_on='pubmed',
                                                     source_processors={'pubmed': [to_list_nan]},
                                                     target_processors={'pubmed': [to_list_nan]})
        source_df = source_df.explode('pubmed')
        target_df = target_df.explode('pubmed')
        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='pubmed', target_on='pubmed')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_connections(df)
