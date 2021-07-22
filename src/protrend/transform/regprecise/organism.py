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

        scientific_names = list(self.genome.loc[:, 'name'])

        dtos = [OrganismDTO(name=name) for name in scientific_names]

        dtos = annotate_organisms(dtos=dtos, names=scientific_names)


