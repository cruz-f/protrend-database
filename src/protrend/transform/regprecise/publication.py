from typing import List

import pandas as pd

from protrend.io import read_json_lines, read_from_stack
from protrend.model import Publication
from protrend.annotation import annotate_publications, PublicationDTO
from protrend.utils.processors import apply_processors, to_int_str
from protrend.transform.regprecise.base import RegPreciseTransformer
from protrend.utils import SetList


class PublicationTransformer(RegPreciseTransformer,
                             source='regprecise',
                             version='0.0.0',
                             node=Publication,
                             order=100,
                             register=True):
    default_transform_stack = {'tf_family': 'TranscriptionFactorFamily.json',
                               'tf': 'TranscriptionFactor.json',
                               'rna': 'RNAFamily.json'}
    columns = SetList(['pmid', 'doi', 'title', 'author', 'year', 'protrend_id'])
    tf_family_columns = SetList(['tffamily_id', 'name', 'url', 'description', 'pubmed', 'regulog'])
    tf_columns = SetList(['collection_id', 'name', 'url', 'description', 'pubmed', 'regulog'])
    rna_columns = SetList(['riboswitch_id', 'name', 'url', 'description', 'pubmed', 'rfam', 'regulog'])

    def _transform_tf_family(self, tf_family: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=tf_family, subset=['name'], perfect_match=True, preserve_nan=False)

        df = df.drop(columns=['tffamily_id', 'name', 'url', 'description', 'regulog'])

        return df

    def _transform_tf(self, tf: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=tf, subset=['name'], perfect_match=True, preserve_nan=False)

        df = df.drop(columns=['collection_id', 'name', 'url', 'description', 'regulog'])
        return df

    def _transform_rna(self, rna: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=rna, subset=['name'], perfect_match=True, preserve_nan=False)

        df = df.drop(columns=['riboswitch_id', 'name', 'url', 'description', 'rfam', 'regulog'])
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
        tf_family = read_from_stack(stack=self.transform_stack, file='tf_family',
                                    default_columns=self.tf_family_columns, reader=read_json_lines)
        tf_family = self._transform_tf_family(tf_family)

        tf = read_from_stack(stack=self.transform_stack, file='tf',
                             default_columns=self.tf_columns, reader=read_json_lines)
        tf = self._transform_tf(tf)

        rna = read_from_stack(stack=self.transform_stack, file='rna',
                              default_columns=self.rna_columns, reader=read_json_lines)
        rna = self._transform_rna(rna)

        df = pd.concat([tf_family, tf, rna])
        df = df.explode('pubmed')
        df = df.dropna(subset=['pubmed'])
        df = self.drop_duplicates(df=df, subset=['pubmed'], preserve_nan=False)

        pmids = df['pubmed'].tolist()
        df = self._transform_publications(pmids)

        df = df.drop(columns=['input_value'])

        df = apply_processors(df, pmid=to_int_str, year=to_int_str)

        self._stack_transformed_nodes(df)

        return df
