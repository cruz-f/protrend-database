from typing import List

import pandas as pd

from protrend.io import read_from_stack, read_txt
from protrend.model.model import Publication
from protrend.transform import PublicationDTO
from protrend.transform.annotation import annotate_publications
from protrend.transform.processors import apply_processors, to_int_str
from protrend.transform.regulondb.base import RegulondbTransformer
from protrend.utils import SetList


class PublicationTransformer(RegulondbTransformer):
    default_node = Publication
    default_transform_stack = {'publication': 'publication.txt'}
    default_order = 100
    columns = SetList(['pmid', 'doi', 'title', 'author', 'year',
                       'publication_id', 'reference_id', 'external_db_id', 'source',
                       'publication_note', 'publication_internal_comment'])
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
        df = df.drop(columns=['author_regulondb', 'year_regulondb', 'title_regulondb'])
        df = df.rename(columns={'author_annotation': 'author',
                                'year_annotation': 'year',
                                'title_regulondb': 'title'})
        df = df.drop(columns=['input_value'])

        df = apply_processors(df, pmid=to_int_str)

        self._stack_transformed_nodes(df)
