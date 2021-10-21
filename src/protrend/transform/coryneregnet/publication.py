from typing import List

import pandas as pd

from protrend.model.model import Publication
from protrend.transform import PublicationDTO
from protrend.transform.annotation import annotate_publications
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer
from protrend.transform.processors import apply_processors, to_int_str
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
                       'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence', 'PMID', 'Source'])

    def _transform_publication(self, regulation: pd.DataFrame) -> pd.DataFrame:
        regulation = self.drop_duplicates(df=regulation, subset=['PMID'], perfect_match=True, preserve_nan=True)
        regulation = regulation.dropna(subset=['PMID'])

        def split_pmid(item: str) -> list:
            items = item.split(',')
            return [item.lstrip().rstrip() for item in items]

        regulation = apply_processors(regulation, PMID=split_pmid)

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