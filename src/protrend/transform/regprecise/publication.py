from typing import List, Set, Union

import pandas as pd

from protrend.io.json import read_json_lines
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

    def _read_tf_family(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('tf_family')

        if file_path:
            df = read_json_lines(file_path)

        else:
            df = pd.DataFrame(columns=['tffamily_id', 'name', 'url', 'description', 'pubmed', 'regulog'])

        return df

    def _read_tf(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('tf')

        if file_path:
            df = read_json_lines(file_path)

        else:
            df = pd.DataFrame(columns=['collection_id', 'name', 'url', 'description', 'pubmed', 'regulog'])

        return df

    def _read_rna(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('rna')

        if file_path:
            df = read_json_lines(file_path)

        else:
            df = pd.DataFrame(columns=['riboswitch_id', 'name', 'url', 'description', 'pubmed', 'rfam', 'regulog'])

        return df

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

        publications = pd.DataFrame([dto.to_dict() for dto in dtos])

        if publications.empty:
            publications = pd.DataFrame(columns=['input_value', 'pmid'])

        apply_processors(nan_to_str, df=publications, col='pmid')

        return publications

    def transform(self):

        tf_family = self._read_tf_family()
        tf_family = self._transform_tf_family(tf_family)

        tf = self._read_tf()
        tf = self._transform_tf(tf)

        rna = self._read_rna()
        rna = self._transform_rna(rna)

        df = pd.concat([tf_family, tf, rna], axis=0)

        identifiers = {identifier for pubmed_row in df['pubmed'] for identifier in pubmed_row}
        publications = self._transform_publications(identifiers)

        df = publications.drop(['input_value'], axis=1)

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df

    def connect(self):
        return
