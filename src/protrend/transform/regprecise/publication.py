from typing import List, Set, Union

import pandas as pd

from protrend.io.utils import read_from_stack
from protrend.transform.annotation.publication import annotate_publications
from protrend.transform.dto import PublicationDTO
from protrend.transform.regprecise.settings import PublicationSettings
from protrend.transform.transformer import DefaultTransformer


class PublicationTransformer(DefaultTransformer):
    default_settings = PublicationSettings
    columns = {'pmid', 'doi', 'title', 'author', 'year'}
    tf_family_columns = {'tffamily_id', 'name', 'url', 'description', 'pubmed', 'regulog'}
    tf_columns = {'collection_id', 'name', 'url', 'description', 'pubmed', 'regulog'}
    rna_columns = {'riboswitch_id', 'name', 'url', 'description', 'pubmed', 'rfam', 'regulog'}

    def _transform_tf_family(self, tf_family: pd.DataFrame) -> pd.DataFrame:

        df = self.drop_duplicates(df=tf_family, subset=['name'], perfect_match=True, preserve_nan=False)

        df = df.drop(columns=['tffamily_id', 'name', 'url', 'description', 'regulog'], axis=1)

        return df

    def _transform_tf(self, tf: pd.DataFrame) -> pd.DataFrame:

        df = self.drop_duplicates(df=tf, subset=['name'], perfect_match=True, preserve_nan=False)

        df = df.drop(columns=['collection_id', 'name', 'url', 'description', 'regulog'], axis=1)
        return df

    def _transform_rna(self, rna: pd.DataFrame) -> pd.DataFrame:

        df = self.drop_duplicates(df=rna, subset=['name'], perfect_match=True, preserve_nan=False)

        df = df.drop(columns=['riboswitch_id', 'name', 'url', 'description', 'rfam' 'regulog'], axis=1)
        return df

    @staticmethod
    def _transform_publications(identifiers: Union[Set[str], List[str]]):

        dtos = [PublicationDTO(input_value=identifier) for identifier in identifiers]
        annotate_publications(dtos=dtos, identifiers=identifiers)

        # pmid: List[str]
        # doi: List[str]
        # title: List[str]
        # author: List[str]
        # year: List[str]
        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):

        tf_family = read_from_stack(tl=self, file='tf_family', json=True, default_columns=self.tf_family_columns)
        tf_family = self._transform_tf_family(tf_family)

        tf = read_from_stack(tl=self, file='tf', json=True, default_columns=self.tf_columns)
        tf = self._transform_tf(tf)

        rna = read_from_stack(tl=self, file='rna', json=True, default_columns=self.rna_columns)
        rna = self._transform_rna(rna)

        df = pd.concat([tf_family, tf, rna], axis=0)

        identifiers = {identifier for pubmed_row in df['pubmed'] for identifier in pubmed_row}
        publications = self._transform_publications(identifiers)

        df = publications.drop(['input_value'], axis=1)

        if df.empty:
            df = self.make_empty_frame()

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df
