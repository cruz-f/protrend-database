from typing import List

import pandas as pd

from protrend.io import read_from_stack, read_json_lines
from protrend.model.model import Publication
from protrend.transform import PublicationDTO
from protrend.transform.annotation import annotate_publications
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.transform.processors import apply_processors, to_int_str
from protrend.utils import SetList


class EvidenceTransformer(DBTBSTransformer):
    default_node = Publication
    default_transform_stack = {'operon': 'Operon.json', 'tfbs': 'TFBS.json'}
    default_order = 100
    columns = SetList(['protrend_id', 'operon', 'tfbs', 'pubmed',
                       'pmid', 'doi', 'title', 'author', 'year'])

    operon_columns = SetList(['name', 'tf', 'url', 'evidence', 'pubmed', 'comment', 'gene', 'tfbs'])
    tfbs_columns = SetList(['identifier', 'url', 'regulation', 'pubmed', 'tf', 'operon', 'gene',
                            'location', 'absolute_position', 'sequence'])

    def _transform_publication(self, operon: pd.DataFrame, tfbs: pd.DataFrame) -> pd.DataFrame:
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
