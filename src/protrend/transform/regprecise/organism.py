from typing import Dict, List, Tuple

import pandas as pd

from protrend.transform.annotation.organism import annotate_organisms
from protrend.transform.dto import OrganismDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, nan_to_str
from protrend.transform.regprecise.settings import OrganismSettings
from protrend.transform.transformer import Transformer


class OrganismTransformer(Transformer):

    def __init__(self, settings: OrganismSettings = None):

        if not settings:
            settings = OrganismSettings()

        super().__init__(settings)

    def read(self, *args, **kwargs):
        self.read_json_lines()

    def _transform_genome(self):
        df = self.get('genome', pd.DataFrame(columns=['name']))

        df = df.drop_duplicates(subset=['name'])

        apply_processors(rstrip, lstrip, df=df, col='name')

        df['input_value'] = df['name']

        return df

    @staticmethod
    def _annotate_organisms(names: List[str]):

        dtos = [OrganismDTO(input_value=name) for name in names]
        annotate_organisms(dtos=dtos, names=names)

        organisms = pd.DataFrame([dto.to_dict() for dto in dtos])

        if organisms.empty:
            organisms = pd.DataFrame(columns=['input_value', 'name'])

        apply_processors(nan_to_str, df=organisms, col='name')

        return organisms

    def transform(self) -> pd.DataFrame:

        genome = self._transform_genome()

        names = tuple(genome['input_value'])

        organisms = self._annotate_organisms(names)

        df = pd.merge(organisms, genome, on='input_value', suffixes=('_annotation', '_regprecise'))

        df['name'] = df['name_annotation'].astype(str) + df['name_regprecise'].astype(str)

        df = df.drop(['input_value', 'name_annotation', 'name_regprecise'], axis=1)

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df
