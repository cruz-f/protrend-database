from typing import List

import pandas as pd

from protrend.io import read_from_stack, read_txt, read_json_frame
from protrend.model.model import Publication, Regulator, TFBS, Promoter, Operon, Gene, Organism
from protrend.transform import PublicationDTO
from protrend.annotation import annotate_publications
from protrend.transform.processors import apply_processors, to_int_str
from protrend.transform.regulondb import OrganismTransformer
from protrend.transform.regulondb.base import RegulondbTransformer, RegulondbConnector
from protrend.transform.regulondb.gene import GeneTransformer
from protrend.transform.regulondb.operon import OperonTransformer
from protrend.transform.regulondb.promoter import PromoterTransformer
from protrend.transform.regulondb.regulator import RegulatorTransformer
from protrend.transform.regulondb.tfbs import TFBSTransformer
from protrend.utils import SetList


class PublicationTransformer(RegulondbTransformer):
    default_node = Publication
    default_transform_stack = {'publication': 'publication.txt'}
    default_order = 100
    columns = SetList(['pmid', 'doi', 'title', 'author', 'year',
                       'publication_id', 'reference_id', 'external_db_id', 'source',
                       'publication_note', 'publication_internal_comment', 'protrend_id'])
    read_columns = SetList(['publication_id', 'reference_id', 'external_db_id', 'author', 'title', 'source',
                            'year', 'publication_note', 'publication_internal_comment'])

    def _transform_publication(self, publication: pd.DataFrame) -> pd.DataFrame:
        publication = apply_processors(df=publication, reference_id=to_int_str)
        publication = publication.dropna(subset=['reference_id'])
        publication = self.drop_duplicates(publication, subset=['reference_id'], perfect_match=True, preserve_nan=True)

        publication = self.create_input_value(publication, col='reference_id')

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
        publication = read_from_stack(stack=self.transform_stack, file='publication', default_columns=self.read_columns,
                                      reader=read_txt, skiprows=36, names=self.read_columns)
        df = self._transform_publication(publication)

        pmids = df['input_value'].tolist()
        publications = self._transform_publications(pmids)

        df = pd.merge(publications, df, on='input_value', suffixes=('_annotation', '_regulondb'))
        df = df.drop(columns=['author_regulondb', 'year_regulondb', 'title_regulondb', 'input_value'])
        df = df.rename(columns={'author_annotation': 'author',
                                'year_annotation': 'year',
                                'title_regulondb': 'title'})

        df = apply_processors(df, pmid=to_int_str)

        self._stack_transformed_nodes(df)

        return df


class PublicationToOrganismConnector(RegulondbConnector):
    default_from_node = Publication
    default_to_node = Organism
    default_connect_stack = {'publication': 'integrated_publication.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)

        organism = read_from_stack(stack=self.connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)

        from_identifiers = publication['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = organism['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToRegulatorConnector(RegulondbConnector):
    default_from_node = Publication
    default_to_node = Regulator
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'regulator': 'integrated_regulator.json',
                             'obj_ev_pub': 'object_ev_method_pub_link.txt'}

    def connect(self):
        publication = read_from_stack(stack=self._connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        publication = publication[['protrend_id', 'publication_id']]
        publication = publication.rename(columns={'protrend_id': 'publication_protrend_id'})

        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = regulator[['protrend_id', 'transcription_factor_id', 'sigma_id', 'srna_id']]
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        obj_ev_pub_cols = ['object_id', 'evidence_id', 'method_id', 'publication_id']
        obj_ev_pub = read_from_stack(stack=self._connect_stack, file='obj_ev_pub',
                                     default_columns=obj_ev_pub_cols, reader=read_txt,
                                     names=obj_ev_pub_cols, skiprows=31)

        regulator_tf = regulator.dropna(subset=['transcription_factor_id'])
        obj_tf = pd.merge(obj_ev_pub, regulator_tf, left_on='object_id', right_on='transcription_factor_id')

        regulator_sigma = regulator.dropna(subset=['sigma_id'])
        obj_sigma = pd.merge(obj_ev_pub, regulator_sigma, left_on='object_id', right_on='sigma_id')

        regulator_srna = regulator.dropna(subset=['srna_id'])
        obj_srna = pd.merge(obj_ev_pub, regulator_srna, left_on='object_id', right_on='srna_id')

        obj_reg = pd.concat([obj_tf, obj_sigma, obj_srna], axis=0)

        reg_pub = pd.merge(obj_reg, publication, on='publication_id')
        reg_pub = reg_pub.drop_duplicates(subset=['publication_protrend_id', 'regulator_protrend_id'])

        from_identifiers = reg_pub['publication_protrend_id'].tolist()
        to_identifiers = reg_pub['regulator_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToTFBSConnector(RegulondbConnector):
    default_from_node = Publication
    default_to_node = TFBS
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'tfbs': 'integrated_tfbs.json',
                             'obj_ev_pub': 'object_ev_method_pub_link.txt'}

    def connect(self):
        publication = read_from_stack(stack=self._connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        publication = publication[['protrend_id', 'publication_id']]
        publication = publication.rename(columns={'protrend_id': 'publication_protrend_id'})

        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = tfbs[['protrend_id', 'site_id']]
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id'})

        obj_ev_pub_cols = ['object_id', 'evidence_id', 'method_id', 'publication_id']
        obj_ev_pub = read_from_stack(stack=self._connect_stack, file='obj_ev_pub',
                                     default_columns=obj_ev_pub_cols, reader=read_txt,
                                     names=obj_ev_pub_cols, skiprows=31)

        obj_site = pd.merge(obj_ev_pub, tfbs, left_on='object_id', right_on='site_id')

        tfbs_pub = pd.merge(obj_site, publication, on='publication_id')
        tfbs_pub = tfbs_pub.drop_duplicates(subset=['publication_protrend_id', 'tfbs_protrend_id'])

        from_identifiers = tfbs_pub['publication_protrend_id'].tolist()
        to_identifiers = tfbs_pub['tfbs_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToPromoterConnector(RegulondbConnector):
    default_from_node = Publication
    default_to_node = Promoter
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'promoter': 'integrated_promoter.json',
                             'obj_ev_pub': 'object_ev_method_pub_link.txt'}

    def connect(self):
        publication = read_from_stack(stack=self._connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        publication = publication[['protrend_id', 'publication_id']]
        publication = publication.rename(columns={'protrend_id': 'publication_protrend_id'})

        promoter = read_from_stack(stack=self._connect_stack, file='promoter',
                                   default_columns=PromoterTransformer.columns, reader=read_json_frame)
        promoter = promoter[['protrend_id', 'promoter_id']]
        promoter = promoter.rename(columns={'protrend_id': 'promoter_protrend_id'})

        obj_ev_pub_cols = ['object_id', 'evidence_id', 'method_id', 'publication_id']
        obj_ev_pub = read_from_stack(stack=self._connect_stack, file='obj_ev_pub',
                                     default_columns=obj_ev_pub_cols, reader=read_txt,
                                     names=obj_ev_pub_cols, skiprows=31)

        obj_promoter = pd.merge(obj_ev_pub, promoter, left_on='object_id', right_on='promoter_id')

        promoter_pub = pd.merge(obj_promoter, publication, on='publication_id')
        promoter_pub = promoter_pub.drop_duplicates(subset=['publication_protrend_id', 'promoter_protrend_id'])

        from_identifiers = promoter_pub['publication_protrend_id'].tolist()
        to_identifiers = promoter_pub['promoter_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToOperonConnector(RegulondbConnector):
    default_from_node = Publication
    default_to_node = Operon
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'operon': 'integrated_operon.json',
                             'obj_ev_pub': 'object_ev_method_pub_link.txt'}

    def connect(self):
        publication = read_from_stack(stack=self._connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        publication = publication[['protrend_id', 'publication_id']]
        publication = publication.rename(columns={'protrend_id': 'publication_protrend_id'})

        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = operon[['protrend_id', 'operon_id']]
        operon = operon.rename(columns={'protrend_id': 'operon_protrend_id'})

        obj_ev_pub_cols = ['object_id', 'evidence_id', 'method_id', 'publication_id']
        obj_ev_pub = read_from_stack(stack=self._connect_stack, file='obj_ev_pub',
                                     default_columns=obj_ev_pub_cols, reader=read_txt,
                                     names=obj_ev_pub_cols, skiprows=31)

        obj_operon = pd.merge(obj_ev_pub, operon, left_on='object_id', right_on='operon_id')

        operon_pub = pd.merge(obj_operon, publication, on='publication_id')
        operon_pub = operon_pub.drop_duplicates(subset=['publication_protrend_id', 'operon_protrend_id'])

        from_identifiers = operon_pub['publication_protrend_id'].tolist()
        to_identifiers = operon_pub['operon_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToGeneConnector(RegulondbConnector):
    default_from_node = Publication
    default_to_node = Gene
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'gene': 'integrated_gene.json',
                             'obj_ev_pub': 'object_ev_method_pub_link.txt'}

    def connect(self):
        publication = read_from_stack(stack=self._connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        publication = publication[['protrend_id', 'publication_id']]
        publication = publication.rename(columns={'protrend_id': 'publication_protrend_id'})

        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = gene[['protrend_id', 'gene_id']]
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id'})

        obj_ev_pub_cols = ['object_id', 'evidence_id', 'method_id', 'publication_id']
        obj_ev_pub = read_from_stack(stack=self._connect_stack, file='obj_ev_pub',
                                     default_columns=obj_ev_pub_cols, reader=read_txt,
                                     names=obj_ev_pub_cols, skiprows=31)

        obj_gene = pd.merge(obj_ev_pub, gene, left_on='object_id', right_on='gene_id')

        gene_pub = pd.merge(obj_gene, publication, on='publication_id')
        gene_pub = gene_pub.drop_duplicates(subset=['publication_protrend_id', 'gene_protrend_id'])

        from_identifiers = gene_pub['publication_protrend_id'].tolist()
        to_identifiers = gene_pub['gene_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
