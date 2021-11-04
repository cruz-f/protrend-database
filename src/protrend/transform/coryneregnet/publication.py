from typing import List

import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Publication, Regulator, Organism, Operon, Gene, TFBS, RegulatoryInteraction
from protrend.transform import PublicationDTO
from protrend.annotation import annotate_publications
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer, CoryneRegNetConnector
from protrend.transform.coryneregnet.operon import OperonTransformer
from protrend.transform.coryneregnet.organism import OrganismTransformer
from protrend.transform.coryneregnet.regulator import RegulatorTransformer
from protrend.transform.coryneregnet.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.transform.processors import apply_processors, to_int_str, to_set_list, to_str
from protrend.utils import SetList


class PublicationTransformer(CoryneRegNetTransformer):
    default_node = Publication
    default_transform_stack = {'bsub': 'bsub_regulation.csv',
                               'cglu': 'cglu_regulation.csv',
                               'ecol': 'ecol_regulation.csv',
                               'mtub': 'mtub_regulation.csv'}
    default_order = 100
    columns = SetList(['protrend_id', 'pmid', 'doi', 'title', 'author', 'year',
                       'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                       'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                       'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence', 'PMID', 'Source', 'taxonomy',
                       'coryneregnet_pmid'])

    def _transform_publication(self, regulation: pd.DataFrame) -> pd.DataFrame:
        regulation = self.drop_duplicates(df=regulation, subset=['PMID'], perfect_match=True, preserve_nan=True)
        regulation = regulation.dropna(subset=['PMID'])

        regulation['coryneregnet_pmid'] = regulation['PMID']

        def split_pmid(item: str) -> list:
            items = item.split(',')
            return [item.lstrip().rstrip() for item in items]

        regulation = apply_processors(regulation, PMID=[to_str, split_pmid])

        regulation = regulation.explode(column='PMID')
        regulation = self.drop_duplicates(df=regulation, subset=['PMID'], perfect_match=True, preserve_nan=True)

        regulation = self.create_input_value(regulation, col='PMID')

        return regulation

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
        regulation = self._build_regulations()

        publication = self._transform_publication(regulation)

        pmids = publication['input_value'].tolist()
        publications = self._transform_publications(pmids)

        df = pd.merge(publications, publication, on='input_value', suffixes=('_annotation', '_coryneregnet'))
        df = df.drop(columns=['input_value'])

        df = apply_processors(df, pmid=to_int_str)

        self._stack_transformed_nodes(df)

        return df


class PublicationToOrganismConnector(CoryneRegNetConnector):
    default_from_node = Publication
    default_to_node = Organism
    default_connect_stack = {'publication': 'integrated_publication.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        publication = apply_processors(publication, taxonomy=to_int_str)

        organism = read_from_stack(stack=self.connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = apply_processors(organism, ncbi_taxonomy=to_int_str)

        df = pd.merge(publication, organism, left_on='taxonomy', right_on='ncbi_taxonomy',
                      suffixes=('_publication', '_organism'))

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['protrend_id_organism'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToRegulatorConnector(CoryneRegNetConnector):
    default_from_node = Publication
    default_to_node = Regulator
    default_connect_stack = {'publication': 'integrated_publication.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)

        regulator = read_from_stack(stack=self.connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)

        df = pd.merge(publication, regulator, on='TF_locusTag', suffixes=('_publication', '_regulator'))

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['protrend_id_regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToOperonConnector(CoryneRegNetConnector):
    default_from_node = Publication
    default_to_node = Operon
    default_connect_stack = {'publication': 'integrated_publication.json', 'operon': 'integrated_operon.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)

        operon = read_from_stack(stack=self.connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = apply_processors(operon, Genes=to_set_list)
        operon = operon.explode(column='Genes')

        df = pd.merge(publication, operon, left_on='TG_locusTag', right_on='Genes',
                      suffixes=('_publication', '_operon'))
        df = PublicationTransformer.drop_duplicates(df=df, subset=['protrend_id_publication', 'protrend_id_operon'],
                                                    perfect_match=True, preserve_nan=True)
        df = df.dropna(subset=['protrend_id_publication', 'protrend_id_operon'])

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['protrend_id_operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToGeneConnector(CoryneRegNetConnector):
    default_from_node = Publication
    default_to_node = Gene
    default_connect_stack = {'publication': 'integrated_publication.json', 'operon': 'integrated_operon.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)

        operon = read_from_stack(stack=self.connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = apply_processors(operon, Genes=to_set_list)
        operon = operon.explode(column='Genes')

        df = pd.merge(publication, operon, left_on='TG_locusTag', right_on='Genes',
                      suffixes=('_publication', '_operon'))
        df = PublicationTransformer.drop_duplicates(df=df, subset=['protrend_id_publication', 'protrend_id_operon'],
                                                    perfect_match=True, preserve_nan=True)
        df = df.dropna(subset=['protrend_id_publication', 'protrend_id_operon'])

        df = df.explode(column='genes')
        df = PublicationTransformer.drop_duplicates(df=df, subset=['protrend_id_publication', 'genes'],
                                                    perfect_match=True, preserve_nan=True)
        df = df.dropna(subset=['protrend_id_publication', 'genes'])

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['genes'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToTFBSConnector(CoryneRegNetConnector):
    default_from_node = Publication
    default_to_node = TFBS
    default_connect_stack = {'publication': 'integrated_publication.json', 'operon': 'integrated_operon.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)

        operon = read_from_stack(stack=self.connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = apply_processors(operon, Genes=to_set_list)
        operon = operon.explode(column='Genes')

        df = pd.merge(publication, operon, left_on='TG_locusTag', right_on='Genes',
                      suffixes=('_publication', '_operon'))
        df = PublicationTransformer.drop_duplicates(df=df, subset=['protrend_id_publication', 'protrend_id_operon'],
                                                    perfect_match=True, preserve_nan=True)
        df = df.dropna(subset=['protrend_id_publication', 'protrend_id_operon'])

        df = df.explode(column='tfbss')
        df = PublicationTransformer.drop_duplicates(df=df, subset=['protrend_id_publication', 'tfbss'],
                                                    perfect_match=True, preserve_nan=True)
        df = df.dropna(subset=['protrend_id_publication', 'tfbss'])

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['tfbss'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToRegulatoryInteractionConnector(CoryneRegNetConnector):
    default_from_node = Publication
    default_to_node = RegulatoryInteraction
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        rin = read_from_stack(stack=self.connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        df = pd.merge(rin, publication, left_on='PMID', right_on='coryneregnet_pmid',
                      suffixes=('_rin', '_publication'))
        df = df.drop_duplicates(subset=['protrend_id_publication', 'protrend_id_rin'])

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['protrend_id_rin'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
