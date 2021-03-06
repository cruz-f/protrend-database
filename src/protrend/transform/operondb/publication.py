import pandas as pd

from protrend.io.utils import read_operon
from protrend.model import Publication, Organism, Gene, Operon
from protrend.transform.mix_ins import PublicationMixIn
from protrend.transform.operondb.base import OperonDBTransformer, OperonDBConnector
from protrend.transform.operondb.operon import OperonTransformer
from protrend.transform.transformations import drop_empty_string, drop_duplicates, create_input_value, merge_columns
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, to_int_str, to_list_nan, to_str


class PublicationTransformer(PublicationMixIn, OperonDBTransformer,
                             source='operondb',
                             version='0.0.0',
                             node=Publication,
                             order=80,
                             register=True):
    columns = SetList(['protrend_id', 'pmid', 'doi', 'title', 'author', 'year',
                       'operon',
                       'operon_db_id', 'name', 'function', 'genes', 'strand', 'start', 'stop',
                       'organism', 'source'])

    @staticmethod
    def transform_publication(operon: pd.DataFrame) -> pd.DataFrame:
        operon = apply_processors(operon, source=to_str)

        # only the known operons have valid pmid
        mask = operon['operon_db_id'].str.startswith('K')
        operon = operon[mask].copy()

        operon = operon.rename(columns={'protrend_id': 'operon'})
        operon = operon.assign(pmid=operon['source'].str.split(' '))

        operon = apply_processors(operon, pmid=to_list_nan)
        operon = operon.explode(column='pmid')

        operon = operon.dropna(subset=['pmid'])
        operon = drop_empty_string(operon, 'pmid')
        operon = drop_duplicates(df=operon, subset=['pmid'])

        operon = create_input_value(operon, col='pmid')
        return operon

    def transform(self):
        operon = read_operon(source=self.source, version=self.version, columns=OperonTransformer.columns)

        # noinspection DuplicatedCode
        publications = self.transform_publication(operon)
        annotated_publications = self.annotate_publications(publications)

        df = pd.merge(annotated_publications, publications, on='input_value', suffixes=('_annotation', '_literature'))
        df = merge_columns(df=df, column='pmid', left='pmid_annotation', right='pmid_literature')

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

    def connect(self):
        df = self.create_connection(source='publication', target='publication', target_column='organism')
        self.stack_connections(df)


class PublicationToOperonConnector(OperonDBConnector,
                                   source='operondb',
                                   version='0.0.0',
                                   from_node=Publication,
                                   to_node=Operon,
                                   register=True):

    def connect(self):
        df = self.create_connection(source='publication', target='publication', target_column='operon')
        self.stack_connections(df)


class PublicationToGeneConnector(OperonDBConnector,
                                 source='operondb',
                                 version='0.0.0',
                                 from_node=Publication,
                                 to_node=Gene,
                                 register=True):

    def connect(self):
        source_df, target_df = self.transform_stacks(source='publication',
                                                     target='gene',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_on='operon_db_id',
                                                     target_on='operon_db_id',
                                                     source_processors={},
                                                     target_processors={'operon_db_id': [to_list_nan]})
        target_df = target_df.explode('operon_db_id')

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='operon_db_id', target_on='operon_db_id')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_connections(df)
