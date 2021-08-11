from typing import List

import pandas as pd

from protrend.io.csv import read_csv
from protrend.io.utils import read_from_stack
from protrend.model.model import Source
from protrend.transform.annotation.pathway import annotate_pathways
from protrend.transform.connector import DefaultConnector
from protrend.transform.dto import PathwayDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, nan_to_str
from protrend.transform.regprecise.settings import PathwaySettings, PathwayToSource
from protrend.transform.regprecise.source import SourceTransformer
from protrend.transform.transformer import DefaultTransformer


class PathwayTransformer(DefaultTransformer):
    default_settings = PathwaySettings
    columns = {'protrend_id',
               'pathway_id', 'name', 'url', 'regulog',
               'synonyms', 'kegg_pathways'}
    read_columns = {'pathway_id', 'name', 'url', 'regulog'}

    def _transform_pathway(self, pathway: pd.DataFrame) -> pd.DataFrame:

        df = self.drop_duplicates(df=pathway, subset=['name'], perfect_match=True, preserve_nan=False)

        apply_processors(rstrip, lstrip, df=df, col='name')

        df['input_value'] = df['name']

        return df

    @staticmethod
    def _transform_pathways(names: List[str]):

        dtos = [PathwayDTO(input_value=name) for name in names]
        annotate_pathways(dtos=dtos, names=names)

        pathways = pd.DataFrame([dto.to_dict() for dto in dtos])

        # name: List[str]
        # synonyms: List[str]
        # kegg_pathways: List[str]

        if pathways.empty:
            pathways = pd.DataFrame(columns=['input_value', 'name'])

        apply_processors(nan_to_str, df=pathways, col='name')

        return pathways

    def transform(self):

        pathway = read_from_stack(tl=self, file='pathway', json=True, default_columns=self.read_columns)
        pathway = self._transform_pathway(pathway)

        names = list(pathway['input_value'])

        pathways = self._transform_pathways(names)

        df = pd.merge(pathways, pathway, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise', fill='')

        df = df.drop(['input_value'], axis=1)

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df


class PathwayToSourceConnector(DefaultConnector):
    default_settings = PathwayToSource

    def connect(self):
        organism = read_from_stack(tl=self, file='pathway', json=False, default_columns=PathwayTransformer.columns)
        source = read_from_stack(tl=self, file='source', json=False, default_columns=SourceTransformer.columns)

        from_identifiers = organism['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=organism['url'].tolist(),
                      external_identifier=organism['pathway_id'].tolist(),
                      key=['pathway_id'] * size)

        df = self.make_connection(size=size,
                                  from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_csv(df)
