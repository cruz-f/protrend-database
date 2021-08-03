from typing import List

import pandas as pd

from protrend.transform.annotation.pathway import annotate_pathways
from protrend.transform.dto import PathwayDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, nan_to_str
from protrend.transform.regprecise.settings import PathwaySettings
from protrend.transform.transformer import Transformer


class PathwayTransformer(Transformer):

    def __init__(self, settings: PathwaySettings = None):

        if not settings:
            settings = PathwaySettings()

        super().__init__(settings)

    def read(self, **kwargs):
        return self._read_json_lines()

    @staticmethod
    def _transform_pathway(df):
        df = df.drop_duplicates(subset=['name'])

        apply_processors(rstrip, lstrip, df=df, col='name')

        df['input_value'] = df['name']

        return df

    @staticmethod
    def _transform_pathways(names: List[str]):

        dtos = [PathwayDTO(input_value=name) for name in names]
        annotate_pathways(dtos=dtos, names=names)

        pathways = pd.DataFrame([dto.to_dict() for dto in dtos])

        if pathways.empty:
            pathways = pd.DataFrame(columns=['input_value', 'name'])

        apply_processors(nan_to_str, df=pathways, col='name')

        return pathways

    def transform(self, **kwargs):

        pathway = kwargs.get('pathway', pd.DataFrame(columns=['name']))
        pathway = self._transform_pathway(pathway)

        names = tuple(pathway['input_value'])

        effectors = self._transform_pathways(names)

        df = pd.merge(effectors, pathway, on='input_value', suffixes=('_annotation', '_regprecise'))

        df['name'] = df['name_annotation'].astype(str) + df['name_regprecise'].astype(str)

        df = df.drop(['input_value', 'name_annotation', 'name_regprecise'], axis=1)

        n_rows, _ = df.shape
        df['source_db'] = ['regprecise'] * n_rows
        df['api_key'] = ['pathway_id'] * n_rows

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df
