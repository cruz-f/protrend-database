from typing import List

import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Publication, Organism, Regulator, Operon, Gene, RegulatoryInteraction
from protrend.transform import PublicationDTO
from protrend.annotation import annotate_publications
from protrend.transform.literature.base import LiteratureTransformer, LiteratureConnector
from protrend.transform.literature.operon import OperonTransformer
from protrend.transform.literature.organism import OrganismTransformer
from protrend.transform.literature.regulator import RegulatorTransformer
from protrend.transform.literature.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.transform.processors import apply_processors, to_int_str, to_set_list
from protrend.utils import SetList


class PublicationTransformer(LiteratureTransformer):
    default_node = Publication
    default_order = 100
    columns = SetList(['protrend_id', 'pmid', 'doi', 'title', 'author', 'year',
                       'regulator_locus_tag', 'operon', 'genes_locus_tag',
                       'regulatory_effect', 'evidence', 'effector', 'mechanism',
                       'publication', 'taxonomy', 'source', 'network_id'])

    def _transform_publication(self, network: pd.DataFrame) -> pd.DataFrame:
        network = apply_processors(network, publication=[to_set_list])
        network = network.explode(column='publication')

        network = self.drop_duplicates(df=network, subset=['publication'], perfect_match=True, preserve_nan=True)
        network = network.dropna(subset=['publication'])

        network = self.create_input_value(network, col='publication')

        return network

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
        network = self._build_network()
        publication = self._transform_publication(network)

        pmids = publication['input_value'].tolist()
        publications = self._transform_publications(pmids)

        df = pd.merge(publications, publication, on='input_value', suffixes=('_annotation', '_literature'))
        df = df.drop(columns=['input_value'])

        df = apply_processors(df, pmid=to_int_str)

        self._stack_transformed_nodes(df)

        return df


class PublicationToOrganismConnector(LiteratureConnector):
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


class PublicationToRegulatorConnector(LiteratureConnector):
    default_from_node = Publication
    default_to_node = Regulator
    default_connect_stack = {'publication': 'integrated_publication.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)

        regulator = read_from_stack(stack=self.connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)

        df = pd.merge(publication, regulator, on='network_id', suffixes=('_publication', '_regulator'))

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['protrend_id_regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToOperonConnector(LiteratureConnector):
    default_from_node = Publication
    default_to_node = Operon
    default_connect_stack = {'publication': 'integrated_publication.json', 'operon': 'integrated_operon.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)

        operon = read_from_stack(stack=self.connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)

        df = pd.merge(publication, operon, on='network_id', suffixes=('_publication', '_operon'))

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['protrend_id_operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToGeneConnector(LiteratureConnector):
    default_from_node = Publication
    default_to_node = Gene
    default_connect_stack = {'publication': 'integrated_publication.json', 'operon': 'integrated_operon.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)

        operon = read_from_stack(stack=self.connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)

        df = pd.merge(publication, operon, on='network_id', suffixes=('_publication', '_operon'))
        df = df.explode(column='genes')

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['genes'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class PublicationToRegulatoryInteractionConnector(LiteratureConnector):
    default_from_node = Publication
    default_to_node = RegulatoryInteraction
    default_connect_stack = {'publication': 'integrated_publication.json',
                             'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        publication = read_from_stack(stack=self.connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        rin = read_from_stack(stack=self.connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        df = pd.merge(publication, rin, on='network_id', suffixes=('_publication', '_rin'))

        from_identifiers = df['protrend_id_publication'].tolist()
        to_identifiers = df['protrend_id_rin'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
