import os
from functools import partial
from typing import Dict

import pandas as pd

from protrend.transform.annotation.organism import annotate_organisms
from protrend.transform.dto import OrganismDTO
from protrend.transform.regprecise.settings import RegPreciseTransformSettings
from protrend.transform.transformer import Transformer


class OrganismTransformer(Transformer):

    def __init__(self,
                 source: str = None,
                 version: str = None,
                 **files: Dict[str, str]):

        if not source:
            source = RegPreciseTransformSettings.source

        if not version:
            version = RegPreciseTransformSettings.version

        if not files:
            files = RegPreciseTransformSettings.organism

        super().__init__(source=source, version=version, **files)

    def read(self, *args, **kwargs):
        self.read_json_lines()

    def process(self):

        self.genome: pd.DataFrame

        scientific_names = self.genome.loc[:, 'name']

        dtos = []
        for name in scientific_names:
            dto = OrganismDTO()
            dto.name.append(name)
            dtos.append(dto)

        dtos = annotate_organisms(dtos=dtos, names=scientific_names)

        dfs = [dto.to_df() for dto in dtos]

        self._df = pd.concat(dfs)

    def integrate(self):
        pass
