from typing import List

import pandas as pd

from protrend.io import read_from_stack, read_json_lines
from protrend.model.model import Publication
from protrend.transform import PublicationDTO
from protrend.transform.annotation import annotate_publications
from protrend.transform.collectf.base import CollectfTransformer
from protrend.transform.processors import apply_processors, parse_collectf_pubmed, to_int_str


class PublicationTransformer(CollectfTransformer):
    default_node = Publication
    default_node_factors = ('pmid', )
    default_transform_stack = {'tf': 'TranscriptionFactor.json'}
    default_order = 100
    columns = {'protrend_id',
               'pmid', 'doi', 'title', 'author', 'year'
               'name', 'family', 'description', 'regulon',
               'pubmed'}
    read_columns = {'name', 'family', 'description', 'regulon'}

    def _transform_tf(self, tf: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=tf, subset=['name'], perfect_match=True, preserve_nan=True)

        df['pubmed'] = df['description']
        df = apply_processors(df=df, pubmed=parse_collectf_pubmed)
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
        tf = read_from_stack(stack=self.transform_stack, file='tf',
                             default_columns=self.read_columns, reader=read_json_lines)

        df = self._transform_tf(tf)
        df = df.dropna(subset=['pubmed'])
        df = self.drop_duplicates(df=df, subset=['pubmed'], perfect_match=False, preserve_nan=False)

        pmids = df['input_value'].tolist()
        publications = self._transform_publications(pmids)

        df = pd.merge(publications, df, on='input_value', suffixes=('_annotation', '_collectf'))

        df = df.drop(columns=['input_value'])

        df = apply_processors(df, pmid=to_int_str)

        self._stack_transformed_nodes(df)

        return df
