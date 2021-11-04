from typing import List

import pandas as pd

from protrend.io import read_from_stack, read_json_lines, read_json_frame
from protrend.model.model import Publication, Regulator, Operon, Gene, TFBS, RegulatoryInteraction
from protrend.transform import PublicationDTO
from protrend.annotation import annotate_publications
from protrend.transform.collectf.base import CollectfTransformer, CollectfConnector
from protrend.transform.collectf.gene import GeneTransformer
from protrend.transform.collectf.operon import OperonTransformer
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.transform.collectf.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.transform.collectf.tfbs import TFBSTransformer
from protrend.utils.processors import apply_processors, to_int_str, to_list
from protrend.utils import SetList


class PublicationTransformer(CollectfTransformer):
    default_node = Publication
    default_transform_stack = {'tfbs': 'TFBS.json'}
    default_order = 100
    columns = SetList(['pmid', 'doi', 'title', 'author', 'year', 'tfbs_id', 'site_start',
                       'site_end', 'site_strand', 'mode', 'sequence', 'pubmed', 'organism',
                       'regulon', 'experimental_evidence', 'operon', 'gene', 'protrend_id'])
    read_columns = SetList(['tfbs_id', 'site_start', 'site_end', 'site_strand', 'mode', 'sequence',
                            'pubmed', 'organism', 'regulon', 'operon', 'gene', 'experimental_evidence'])

    def _transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        df = apply_processors(df=tfbs, pubmed=to_list)
        df = df.explode(column='pubmed')

        df = self.create_input_value(df, col='pubmed')

        return df

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
        tfbs = read_from_stack(stack=self.transform_stack, file='tfbs',
                               default_columns=self.read_columns, reader=read_json_lines)

        df = self._transform_tfbs(tfbs)
        df = df.dropna(subset=['pubmed'])
        df = self.drop_duplicates(df=df, subset=['pubmed'], perfect_match=False, preserve_nan=False)

        pmids = df['input_value'].tolist()
        publications = self._transform_publications(pmids)

        df = pd.merge(publications, df, on='input_value', suffixes=('_annotation', '_collectf'))

        df = df.drop(columns=['input_value'])

        df = apply_processors(df, pmid=to_int_str)

        self._stack_transformed_nodes(df)

        return df


class PublicationToRegulatorConnector(CollectfConnector):
    default_from_node = Publication
    default_to_node = Regulator
    default_connect_stack = {'publication': 'integrated_publication.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        publication = apply_processors(publication, regulon=to_list)
        publication = publication.explode(column='regulon')
        publication = publication.rename(columns={'protrend_id': 'publication_protrend_id'})

        regulator = read_from_stack(stack=self.connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        df = pd.merge(publication, regulator, left_on='regulon', right_on='uniprot_accession')
        df = df.dropna(subset=['publication_protrend_id', 'regulator_protrend_id'])
        df = df.drop_duplicates(subset=['publication_protrend_id', 'regulator_protrend_id'])

        from_identifiers = df['publication_protrend_id'].tolist()
        to_identifiers = df['regulator_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)
        self.stack_json(df)


class PublicationToOperonConnector(CollectfConnector):
    default_from_node = Publication
    default_to_node = Operon
    default_connect_stack = {'publication': 'integrated_publication.json', 'operon': 'integrated_operon.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        publication = apply_processors(publication, operon=to_list)
        publication = publication.explode(column='operon')
        publication = publication.rename(columns={'protrend_id': 'publication_protrend_id'})

        operon = read_from_stack(stack=self.connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = apply_processors(operon, operon_id_old=to_list)
        operon = operon.explode(column='operon_id_old')
        operon = operon.rename(columns={'protrend_id': 'operon_protrend_id'})

        df = pd.merge(publication, operon, left_on='operon', right_on='operon_id_old')
        df = df.dropna(subset=['publication_protrend_id', 'operon_protrend_id'])
        df = df.drop_duplicates(subset=['publication_protrend_id', 'operon_protrend_id'])

        from_identifiers = df['publication_protrend_id'].tolist()
        to_identifiers = df['operon_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)
        self.stack_json(df)


class PublicationToGeneConnector(CollectfConnector):
    default_from_node = Publication
    default_to_node = Gene
    default_connect_stack = {'publication': 'integrated_publication.json', 'gene': 'integrated_gene.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        publication = apply_processors(publication, gene=to_list)
        publication = publication.explode(column='gene')
        publication = publication.rename(columns={'protrend_id': 'publication_protrend_id'})

        gene = read_from_stack(stack=self.connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = apply_processors(gene, locus_tag_old=to_list)
        gene = gene.explode(column='locus_tag_old')
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id'})

        df = pd.merge(publication, gene, left_on='gene', right_on='locus_tag_old')
        df = df.dropna(subset=['publication_protrend_id', 'gene_protrend_id'])
        df = df.drop_duplicates(subset=['publication_protrend_id', 'gene_protrend_id'])

        from_identifiers = df['publication_protrend_id'].tolist()
        to_identifiers = df['gene_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)
        self.stack_json(df)


class PublicationToTFBSConnector(CollectfConnector):
    default_from_node = Publication
    default_to_node = TFBS
    default_connect_stack = {'publication': 'integrated_publication.json', 'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        publication = publication.rename(columns={'protrend_id': 'publication_protrend_id'})

        tfbs = read_from_stack(stack=self.connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id'})

        df = pd.merge(publication, tfbs, on='tfbs_id')
        df = df.dropna(subset=['publication_protrend_id', 'tfbs_protrend_id'])
        df = df.drop_duplicates(subset=['publication_protrend_id', 'tfbs_protrend_id'])

        from_identifiers = df['publication_protrend_id'].tolist()
        to_identifiers = df['tfbs_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)
        self.stack_json(df)


class PublicationToRegulatoryInteractionConnector(CollectfConnector):
    default_from_node = Publication
    default_to_node = RegulatoryInteraction
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        publication = apply_processors(publication, regulon=to_list)
        publication = publication.explode(column='regulon')
        publication = publication.rename(columns={'protrend_id': 'publication_protrend_id'})

        rin = read_from_stack(stack=self.connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        rin = apply_processors(rin, regulon=to_list)
        rin = rin.explode(column='regulon')
        rin = rin.rename(columns={'protrend_id': 'rin_protrend_id'})

        df = pd.merge(publication, rin, on='regulon')
        df = df.dropna(subset=['publication_protrend_id', 'rin_protrend_id'])
        df = df.drop_duplicates(subset=['publication_protrend_id', 'rin_protrend_id'])

        from_identifiers = df['publication_protrend_id'].tolist()
        to_identifiers = df['rin_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)
        self.stack_json(df)
