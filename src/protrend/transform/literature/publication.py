from typing import List

import pandas as pd

from protrend.model.model import Publication
from protrend.transform import PublicationDTO
from protrend.transform.annotation import annotate_publications
from protrend.transform.literature.base import LiteratureTransformer
from protrend.transform.processors import apply_processors, to_int_str, to_set_list
from protrend.utils import SetList


class PublicationTransformer(LiteratureTransformer):
    default_node = Publication
    default_order = 100
    columns = SetList(['protrend_id', 'pmid', 'doi', 'title', 'author', 'year',
                       'regulator_locus_tag', 'regulator_name', 'operon', 'genes_locus_tag',
                       'genes_name', 'regulatory_effect', 'evidence', 'effector', 'mechanism',
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
