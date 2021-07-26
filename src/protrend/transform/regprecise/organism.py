import os
from functools import partial
from typing import Dict

import numpy as np
import pandas as pd

from protrend.model.model import Organism
from protrend.transform.annotation.organism import annotate_organisms
from protrend.transform.dto import OrganismDTO
from protrend.transform.regprecise.settings import RegPreciseTransformSettings
from protrend.transform.transformer import Transformer


class OrganismTransformer(Transformer):

    node = Organism

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

    @property
    def df(self) -> pd.DataFrame:

        if self._df.empty:
            return pd.DataFrame(columns=list(self.node.cls_keys()))

        return self._df

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
        snapshot = self.node_snapshot()

        create = []
        update = []
        remove = []

        for i, organism in self._df.iterrows():

            if organism['ncbi_taxonomy']:

                if organism['ncbi_taxonomy'] in snapshot.loc['ncbi_taxonomy']:
                    update.append(i)

                else:
                    create.append(i)

            elif organism['name']:

                if organism['name'] in snapshot.loc['name']:
                    update.append(i)

                else:
                    create.append(i)

            else:
                remove.append(i)

