import pandas as pd

from protrend import Operon
from protrend.io import read_from_stack, read_json_frame
from protrend.model import Publication, Organism, Gene
from protrend.transform.operondb.base import OperonDBTransformer, OperonDBConnector
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, to_int_str, to_list_nan


class PublicationTransformer(OperonDBTransformer,
                             source='operondb',
                             version='0.0.0',
                             node=Publication,
                             order=80,
                             register=True):
    default_transform_stack = {'operon': 'integrated_operon.json'}
    columns = SetList(['protrend_id', 'pmid', 'doi', 'title', 'author', 'year',
                       'operon',
                       'operon_db_id', 'name', 'function', 'genes', 'strand', 'start', 'stop',
                       'organism', 'pubmed'])

    def transform_publication(self, operon: pd.DataFrame) -> pd.DataFrame:
        operon = operon.rename(columns={'protrend_id': 'operon'})
        operon = operon.assign(pmid=operon['pubmed'].str.split(' '))

        operon = apply_processors(operon, pmid=to_list_nan)
        operon = operon.explode(column='pmid')

        operon = operon.dropna(subset=['pmid'])
        operon = self.drop_empty_string(operon, 'pmid')
        operon = self.drop_duplicates(df=operon, subset=['pmid'])

        operon = self.create_input_value(operon, col='pmid')
        return operon

    def transform(self):
        operon = read_from_stack(stack=self.transform_stack, key='operon', columns=self.read_columns,
                                 reader=read_json_frame)

        # noinspection DuplicatedCode
        publications = self.transform_publication(operon)
        annotated_publications = self.annotate_publications(publications)

        df = pd.merge(annotated_publications, publications, on='input_value', suffixes=('_annotation', '_literature'))
        df = apply_processors(df, pmid=to_int_str, year=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df


class PublicationToOrganismConnector(OperonDBConnector,
                                     source='operondb',
                                     version='0.0.0',
                                     from_node=Publication,
                                     to_node=Organism,
                                     register=True):
    default_connect_stack = {'publication': 'integrated_publication.json'}

    def connect(self):
        df = self.create_connection(source='publication', target='publication', target_column='organism')
        self.stack_json(df)


class PublicationToOperonConnector(OperonDBConnector,
                                   source='operondb',
                                   version='0.0.0',
                                   from_node=Publication,
                                   to_node=Operon,
                                   register=True):
    default_connect_stack = {'publication': 'integrated_publication.json'}

    def connect(self):
        df = self.create_connection(source='publication', target='publication', target_column='operon')
        self.stack_json(df)


class PublicationToGeneConnector(OperonDBConnector,
                                 source='operondb',
                                 version='0.0.0',
                                 from_node=Publication,
                                 to_node=Gene,
                                 register=True):
    default_connect_stack = {'publication': 'integrated_publication.json'}

    def connect(self):
        source_df, target_df = self.transform_stacks(source='publication',
                                                     target='publication',
                                                     source_column='protrend_id',
                                                     target_column='genes',
                                                     source_processors={},
                                                     target_processors={'genes': [to_list_nan]})
        target_df = target_df.explode('genes')

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='operon_db_id', target_on='operon_db_id')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_json(df)
