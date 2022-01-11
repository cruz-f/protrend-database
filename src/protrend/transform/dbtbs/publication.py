import pandas as pd

from protrend.io import read_from_stack, read_json_lines
from protrend.model import Organism, Publication, Regulator, Gene, TFBS, RegulatoryInteraction
from protrend.transform.dbtbs.base import DBTBSTransformer, DBTBSConnector
from protrend.transform.dbtbs.tfbs import TFBSTransformer
from protrend.transform.mix_ins import PublicationMixIn
from protrend.transform.transformations import (drop_empty_string, select_columns, group_by, create_input_value,
                                                merge_columns)
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, to_int_str, to_list_nan, to_set_list, flatten_set_list


class PublicationTransformer(PublicationMixIn, DBTBSTransformer,
                             source='dbtbs',
                             version='0.0.4',
                             node=Publication,
                             order=100,
                             register=True):
    default_transform_stack = {'tfbs': 'TFBS.json'}
    columns = SetList(['protrend_id', 'pmid', 'doi', 'title', 'author', 'year',
                       'pubmed', 'tf', 'gene', 'tfbs'])

    @staticmethod
    def transform_publication(tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = tfbs.assign(pmid=tfbs['pubmed'].copy(), tfbs=tfbs['identifier'].copy())
        tfbs = apply_processors(tfbs, pmid=to_list_nan)

        tfbs = tfbs.explode(column='pmid')
        tfbs = apply_processors(tfbs, pmid=to_int_str)

        tfbs = tfbs.dropna(subset=['pmid'])
        tfbs = drop_empty_string(tfbs, 'pmid')
        tfbs = select_columns(tfbs, 'pubmed', 'pmid', 'tf', 'gene', 'tfbs')

        aggregation = {'tfbs': to_set_list}
        tfbs = group_by(tfbs, column='pmid', aggregation=aggregation, default=flatten_set_list)

        tfbs = create_input_value(tfbs, col='pmid')
        return tfbs

    def transform(self):
        tfbs = read_from_stack(stack=self.transform_stack, key='tfbs',
                               columns=TFBSTransformer.columns, reader=read_json_lines)

        publications = self.transform_publication(tfbs)
        annotated_publications = self.annotate_publications(publications)

        df = pd.merge(annotated_publications, publications, on='input_value', suffixes=('_annotation', '_dbtbs'))

        df = merge_columns(df=df, column='pmid', left='pmid_annotation', right='pmid_dbtbs')
        df = apply_processors(df, pmid=to_int_str, year=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df


class PublicationToOrganismConnector(DBTBSConnector,
                                     source='dbtbs',
                                     version='0.0.4',
                                     from_node=Publication,
                                     to_node=Organism,
                                     register=True):
    default_connect_stack = {'publication': 'integrated_publication.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        df = self.create_connection(source='publication', target='organism',
                                    cardinality='many_to_one')
        self.stack_json(df)


class PublicationConnector(DBTBSConnector, register=False):
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'rin': 'integrated_regulatoryinteraction.json'}

    def _connect(self, target_column: str, source_on: str, target_on: str):
        source_df, target_df = self.transform_stacks(source='publication',
                                                     target='rin',
                                                     source_column='protrend_id',
                                                     target_column=target_column,
                                                     source_processors={source_on: [to_list_nan]},
                                                     target_processors={})
        source_df = source_df.explode(source_on)

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on=source_on, target_on=target_on)

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_json(df)


class PublicationToRegulatorConnector(PublicationConnector,
                                      source='dbtbs',
                                      version='0.0.4',
                                      from_node=Publication,
                                      to_node=Regulator,
                                      register=True):

    def connect(self):
        return self._connect(target_column='regulator', source_on='tf', target_on='tf')


class PublicationToGeneConnector(PublicationConnector,
                                 source='dbtbs',
                                 version='0.0.4',
                                 from_node=Publication,
                                 to_node=Gene,
                                 register=True):

    def connect(self):
        return self._connect(target_column='gene', source_on='gene', target_on='tg_gene')


class PublicationToTFBSConnector(PublicationConnector,
                                 source='dbtbs',
                                 version='0.0.4',
                                 from_node=Publication,
                                 to_node=TFBS,
                                 register=True):
    def connect(self):
        return self._connect(target_column='tfbs', source_on='tfbs', target_on='identifier')


class PublicationToRegulatoryInteractionConnector(PublicationConnector,
                                                  source='dbtbs',
                                                  version='0.0.4',
                                                  from_node=Publication,
                                                  to_node=RegulatoryInteraction,
                                                  register=True):
    def connect(self):
        return self._connect(target_column='protrend_id', source_on='tfbs', target_on='identifier')
