from typing import Union

import pandas as pd

from protrend.io import read
from protrend.model import Publication, Regulator, TFBS, Gene, Organism, RegulatoryInteraction
from protrend.transform.mix_ins import PublicationMixIn
from protrend.transform.regulondb.base import RegulonDBTransformer, RegulonDBConnector, regulondb_reader
from protrend.transform.transformations import (select_columns, drop_empty_string, drop_duplicates, create_input_value,
                                                merge_columns)
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, to_int_str


def _validate_pmid(item: str) -> Union[str, None]:
    try:
        return str(int(item))
    except ValueError:
        return


class PublicationTransformer(PublicationMixIn, RegulonDBTransformer,
                             source='regulondb',
                             version='0.0.0',
                             node=Publication,
                             order=100,
                             register=True):
    columns = SetList(['protrend_id', 'pmid', 'doi', 'title', 'author', 'year',
                       'publication_id', 'reference_id', 'external_db_id', 'source'])

    @staticmethod
    def transform_publication(publication: pd.DataFrame) -> pd.DataFrame:
        publication = select_columns(publication, 'publication_id', 'reference_id', 'external_db_id', 'source')
        publication = publication.assign(pmid=publication['reference_id'].copy())

        publication = apply_processors(df=publication, pmid=to_int_str)

        publication = publication.dropna(subset=['pmid'])
        publication = drop_empty_string(publication, 'pmid')
        publication = drop_duplicates(publication, subset=['pmid'])

        publication = apply_processors(publication, pmid=_validate_pmid)
        publication = publication.dropna(subset=['pmid'])

        publication = create_input_value(publication, col='pmid')
        return publication

    def transform(self):
        columns = ['publication_id', 'reference_id', 'external_db_id', 'author', 'title', 'source',
                   'year', 'publication_note', 'publication_internal_comment']
        reader = regulondb_reader(skiprows=36, names=columns)
        publication = read(source=self.source, version=self.version,
                           file='publication.txt', reader=reader,
                           default=pd.DataFrame(columns=columns))

        publications = self.transform_publication(publication)
        annotated_publications = self.annotate_publications(publications)

        df = pd.merge(annotated_publications, publications, on='input_value', suffixes=('_annotation', '_regulondb'))

        # merge pmid
        df = merge_columns(df=df, column='pmid', left='pmid_annotation', right='pmid_regulondb')

        df = apply_processors(df, pmid=to_int_str)

        self.stack_transformed_nodes(df)
        return df


class PublicationToOrganismConnector(RegulonDBConnector,
                                     source='regulondb',
                                     version='0.0.0',
                                     from_node=Publication,
                                     to_node=Organism,
                                     register=True):

    def connect(self):
        df = self.create_connection(source='publication', target='organism', cardinality='many_to_one')
        self.stack_connections(df)


class PublicationConnector(RegulonDBConnector, register=False):

    def _connect(self, target: str, obj_ev_pub_col: str, target_col: str):
        source_df, target_df = self.transform_stacks(source='publication',
                                                     target=target,
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_on='publication_id',
                                                     target_on='publication_id',
                                                     source_processors={},
                                                     target_processors={target_col: []})

        obj_ev_pub_cols = ['object_id', 'evidence_id', 'method_id', 'publication_id']
        obj_ev_pub_reader = regulondb_reader(skiprows=31, names=obj_ev_pub_cols)
        obj_ev_pub = read(source=self.source, version=self.version, file='object_ev_method_pub_link.txt',
                          reader=obj_ev_pub_reader, default=pd.DataFrame(columns=obj_ev_pub_cols))

        target_df = pd.merge(obj_ev_pub, target_df, left_on=obj_ev_pub_col, right_on=target_col)

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='publication_id', target_on='publication_id')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        return df


class PublicationToRegulatorConnector(PublicationConnector,
                                      source='regulondb',
                                      version='0.0.0',
                                      from_node=Publication,
                                      to_node=Regulator,
                                      register=True):

    def connect(self):
        df = self._connect(target='regulator', obj_ev_pub_col='object_id', target_col='regulator_id')
        self.stack_connections(df)


class PublicationToGeneConnector(PublicationConnector,
                                 source='regulondb',
                                 version='0.0.0',
                                 from_node=Publication,
                                 to_node=Gene,
                                 register=True):

    def connect(self):
        df = self._connect(target='gene', obj_ev_pub_col='object_id', target_col='gene_id')
        self.stack_connections(df)


class PublicationToTFBSConnector(PublicationConnector,
                                 source='regulondb',
                                 version='0.0.0',
                                 from_node=Publication,
                                 to_node=TFBS,
                                 register=True):

    def connect(self):
        df = self._connect(target='tfbs', obj_ev_pub_col='object_id', target_col='site_id')
        self.stack_connections(df)


class PublicationToRegulatoryInteractionConnector(PublicationConnector,
                                                  source='regulondb',
                                                  version='0.0.0',
                                                  from_node=Publication,
                                                  to_node=RegulatoryInteraction,
                                                  register=True):

    def connect(self):
        regulator_df = self._connect(target='regulator', obj_ev_pub_col='object_id', target_col='regulator_id')
        gene_df = self._connect(target='gene', obj_ev_pub_col='object_id', target_col='gene_id')
        tfbs_df = self._connect(target='tfbs', obj_ev_pub_col='object_id', target_col='site_id')
        df = pd.concat([regulator_df, gene_df, tfbs_df])
        df = df.reset_index(drop=True)
        self.stack_connections(df)
