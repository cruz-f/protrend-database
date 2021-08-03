from typing import List

import pandas as pd

from protrend.transform.annotation.publication import annotate_publications
from protrend.transform.dto import PublicationDTO
from protrend.transform.processors import remove_white_space, remove_regprecise_more, remove_multiple_white_space, \
    rstrip, lstrip, remove_pubmed, apply_processors, nan_to_str
from protrend.transform.regprecise.settings import PublicationSettings
from protrend.transform.transformer import Transformer


class PublicationTransformer(Transformer):

    def __init__(self, settings: PublicationSettings = None):

        if not settings:
            settings = PublicationSettings()

        super().__init__(settings)

    def read(self, **kwargs):
        return self._read_json_lines()

    @staticmethod
    def _transform_tf_family(df):
        df = df.drop_duplicates(subset=['name'])

        df = df.drop(['description'
                      'regulog'
                      'url'],
                     axis=1)

        apply_processors(remove_white_space,
                         df=df,
                         col='name')

        apply_processors(remove_regprecise_more,
                         remove_pubmed,
                         remove_multiple_white_space,
                         rstrip,
                         lstrip,
                         df=df,
                         col='description')

        return df

    @staticmethod
    def _transform_tf(df):

        df = df.drop_duplicates(subset=['name'])

        apply_processors(remove_white_space,
                         df=df,
                         col='name')

        apply_processors(remove_regprecise_more,
                         remove_pubmed,
                         remove_multiple_white_space,
                         rstrip,
                         lstrip,
                         df=df,
                         col='description')

        return df

    @staticmethod
    def _transform_rna(df):
        df = df.drop_duplicates(subset=['name'])

        apply_processors(remove_white_space,
                         df=df,
                         col='name')

        apply_processors(remove_regprecise_more,
                         remove_pubmed,
                         remove_multiple_white_space,
                         rstrip,
                         lstrip,
                         df=df,
                         col='description')

        return df

    @staticmethod
    def _transform_publications(identifiers: List[str]):

        dtos = [PublicationDTO(input_value=identifier) for identifier in identifiers]
        annotate_publications(dtos=dtos, identifiers=identifiers)

        publications = pd.DataFrame([dto.to_dict() for dto in dtos])

        if publications.empty:
            publications = pd.DataFrame(columns=['input_value', 'pmid'])

        apply_processors(nan_to_str, df=publications, col='pmid')

        return publications

    def transform(self, **kwargs):

        tf_family = kwargs.get('tf_family', pd.DataFrame(columns=['name']))
        tf_family = self._transform_tf_family(tf_family)

        tf = kwargs.get('tf', pd.DataFrame(columns=['name']))
        tf = self._transform_tf(tf)

        rna = kwargs.get('rna', pd.DataFrame(columns=['name']))
        rna = self._transform_rna(rna)

        df = pd.merge(tf_family, tf, on='name', suffixes=('_tf_family', '_tf'))

        # concat pubmed
        df['pubmed'] = df['pubmed_tf_family'].astype(set) + df['pubmed_tf'].astype(set)

        df = df.drop(['description_tf_family', 'description_tf',
                      'regulog_tf_family', 'regulog_tf',
                      'pubmed_tf_family', 'pubmed_tf',
                      'url_tf_family', 'url_tf'], axis=1)

        df = pd.merge(df, rna, on='name', suffixes=('_tf', '_rna'))

        # concat pubmed
        df['pubmed'] = df['pubmed_tf'].astype(set) + df['pubmed_rna'].astype(set)

        df = df.drop(['description_tf', 'description_rna',
                      'regulog_tf', 'regulog_rna',
                      'pubmed_tf', 'pubmed_rna',
                      'url'],
                     axis=1)

        identifiers = [identifier for pubmed_row in df['pubmed'] for identifier in pubmed_row]
        publications = self._transform_publications(identifiers)

        merged_df = pd.merge(df, rna, on='name', suffixes=('_tf', '_rna'))

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df
