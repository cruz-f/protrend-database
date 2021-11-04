from typing import List

import pandas as pd

from protrend.io import read_from_stack, read_json_lines, read_json_frame
from protrend.model.model import Publication, Regulator, Operon, Gene, TFBS, RegulatoryInteraction
from protrend.transform import PublicationDTO
from protrend.annotation import annotate_publications
from protrend.transform.dbtbs.base import DBTBSTransformer, DBTBSConnector
from protrend.transform.dbtbs.operon import OperonTransformer
from protrend.transform.dbtbs.regulator import RegulatorTransformer
from protrend.transform.dbtbs.tfbs import TFBSTransformer
from protrend.transform.processors import apply_processors, to_int_str, to_list
from protrend.utils import SetList


class PublicationTransformer(DBTBSTransformer):
    default_node = Publication
    default_transform_stack = {'operon': 'Operon.json', 'tfbs': 'TFBS.json'}
    default_order = 100
    columns = SetList(['protrend_id', 'operon', 'tfbs', 'pubmed',
                       'pmid', 'doi', 'title', 'author', 'year'])

    operon_columns = SetList(['name', 'tf', 'url', 'evidence', 'pubmed', 'comment', 'gene', 'tfbs'])
    tfbs_columns = SetList(['identifier', 'url', 'regulation', 'pubmed', 'tf', 'operon', 'gene',
                            'location', 'absolute_position', 'sequence'])

    def _transform_publication(self, operon: pd.DataFrame, tfbs: pd.DataFrame) -> pd.DataFrame:
        operon = operon.explode(column='name')
        operon = operon.dropna(subset=['pubmed'])
        operon = operon.explode(column='pubmed')
        operon = self.drop_duplicates(df=operon, subset=['pubmed'], perfect_match=True, preserve_nan=True)
        operon = operon.rename(columns={'name': 'operon'})
        operon = self.select_columns(operon, 'operon', 'pubmed')

        tfbs = tfbs.dropna(subset=['pubmed'])
        tfbs = tfbs.explode(column='pubmed')
        tfbs = self.drop_duplicates(df=tfbs, subset=['pubmed'], perfect_match=True, preserve_nan=True)
        tfbs = tfbs.rename(columns={'identifier': 'tfbs'})
        tfbs = self.select_columns(tfbs, 'tfbs', 'pubmed')

        publication = pd.concat([operon, tfbs], axis=0)
        publication = apply_processors(df=publication, pubmed=to_int_str)

        publication = self.drop_duplicates(df=publication, subset=['pubmed'], perfect_match=True, preserve_nan=True)
        publication = self.create_input_value(publication, col='pubmed')

        return publication

    @staticmethod
    def _transform_publications(identifiers: List[str]):
        dtos = [PublicationDTO(input_value=identifier) for identifier in identifiers]
        annotate_publications(dtos=dtos, identifiers=identifiers)

        # pmid: List[str]
        # doi: List[str]
        # title: List[str]
        # author: List[str]
        # year: List[str]
        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=self.operon_columns, reader=read_json_lines)

        tfbs = read_from_stack(stack=self.transform_stack, file='tfbs',
                               default_columns=self.tfbs_columns, reader=read_json_lines)

        publication = self._transform_publication(operon, tfbs)

        pmids = publication['input_value'].tolist()
        publications = self._transform_publications(pmids)

        df = pd.merge(publications, publication, on='input_value', suffixes=('_annotation', '_dbtbs'))
        df = df.drop(columns=['input_value'])
        df = apply_processors(df, pmid=to_int_str)

        self._stack_transformed_nodes(df)
        return df


class PublicationToRegulatorConnector(DBTBSConnector):
    default_from_node = Publication
    default_to_node = Regulator
    default_connect_stack = {'publication': 'integrated_publication.json', 'regulator': 'integrated_regulator.json',
                             'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        publication = apply_processors(publication, pmid=to_int_str)

        regulator = read_from_stack(stack=self.connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)

        tfbs = read_from_stack(stack=self.connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = tfbs.dropna(subset=['pubmed', 'tf'])
        tfbs = apply_processors(tfbs, pubmed=to_list, tf=to_list)
        tfbs = tfbs.explode(column='pubmed')
        tfbs = tfbs.explode(column='tf')
        tfbs = apply_processors(tfbs, pubmed=to_int_str)

        tfbs_pub = pd.merge(tfbs, publication, left_on='pubmed', right_on='pmid', suffixes=('_tfbs', '_publication'))

        reg_pub = pd.merge(tfbs_pub, regulator, left_on='tf', right_on='name_dbtbs')
        reg_pub = reg_pub.drop_duplicates(subset=['protrend_id_publication', 'protrend_id'])

        from_identifiers = reg_pub['protrend_id_publication'].tolist()
        to_identifiers = reg_pub['protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToOperonConnector(DBTBSConnector):
    default_from_node = Publication
    default_to_node = Operon
    default_connect_stack = {'publication': 'integrated_publication.json', 'operon': 'integrated_operon.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        operon = read_from_stack(stack=self.connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)

        df = pd.merge(operon, publication, left_on='name', right_on='operon', suffixes=('_operon', '_publication'))
        df = df.drop_duplicates(subset=['protrend_id_publication', 'protrend_id_operon'])

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['protrend_id_operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToGeneConnector(DBTBSConnector):
    default_from_node = Publication
    default_to_node = Gene
    default_connect_stack = {'publication': 'integrated_publication.json', 'operon': 'integrated_operon.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        operon = read_from_stack(stack=self.connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)

        df = pd.merge(operon, publication, left_on='name', right_on='operon', suffixes=('_operon', '_publication'))
        df = df.drop_duplicates(subset=['protrend_id_publication', 'protrend_id_operon'])

        df = df.explode(column='genes')
        df = df.drop_duplicates(subset=['protrend_id_publication', 'genes'])

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['genes'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToTFBSConnector(DBTBSConnector):
    default_from_node = Publication
    default_to_node = TFBS
    default_connect_stack = {'publication': 'integrated_publication.json', 'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        tfbs = read_from_stack(stack=self.connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)

        df = pd.merge(tfbs, publication, left_on='identifier', right_on='tfbs', suffixes=('_tfbs', '_publication'))
        df = df.drop_duplicates(subset=['protrend_id_publication', 'protrend_id_tfbs'])

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['protrend_id_tfbs'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToRegulatoryInteractionConnector(DBTBSConnector):
    default_from_node = Publication
    default_to_node = RegulatoryInteraction
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        rin = read_from_stack(stack=self.connect_stack, file='rin',
                              default_columns=TFBSTransformer.columns, reader=read_json_frame)

        df = pd.merge(rin, publication, left_on='identifier', right_on='tfbs', suffixes=('_rin', '_publication'))
        df = df.drop_duplicates(subset=['protrend_id_publication', 'protrend_id_rin'])

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['protrend_id_rin'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
