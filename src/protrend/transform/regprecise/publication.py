from typing import List, Set, Union

import pandas as pd

from protrend.transform.annotation.publication import annotate_publications
from protrend.transform.dto import PublicationDTO
from protrend.transform.processors import apply_processors, nan_to_str
from protrend.transform.regprecise.settings import PublicationSettings
from protrend.transform.transformer import Transformer


class PublicationTransformer(Transformer):

    def __init__(self, settings: PublicationSettings = None):

        if not settings:
            settings = PublicationSettings()

        super().__init__(settings)

    def read(self, **kwargs):
        files = self._read_json_lines()

        return super(PublicationTransformer, self).read(**files)

    @staticmethod
    def _transform_tf_family(df):
        df = df.drop_duplicates(subset=['name'])

        df = df.drop(['tffamily_id',
                      'name',
                      'url',
                      'description',
                      'regulog'],
                     axis=1)

        return df

    @staticmethod
    def _transform_tf(df):
        df = df.drop_duplicates(subset=['name'])

        df = df.drop(['collection_id',
                      'name',
                      'url',
                      'description',
                      'regulog'],
                     axis=1)
        return df

    @staticmethod
    def _transform_rna(df):
        df = df.drop_duplicates(subset=['name'])

        df = df.drop(['riboswitch_id',
                      'name',
                      'url',
                      'description',
                      'rfam'
                      'regulog'],
                     axis=1)
        return df

    @staticmethod
    def _transform_publications(identifiers: Union[Set[str], List[str]]):

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

        df = pd.concat([tf_family, tf, rna], axis=0)

        identifiers = {identifier for pubmed_row in df['pubmed'] for identifier in pubmed_row}
        publications = self._transform_publications(identifiers)

        df = publications.drop(['input_value'], axis=1)

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df

    def connect(self, df: pd.DataFrame):
        return
