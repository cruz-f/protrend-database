from typing import Dict

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
        self.taxonomy: pd.DataFrame

        scientific_names = self.genome.loc[:, 'name']

        dtos = []
        for name in scientific_names:
            dto = OrganismDTO()
            dto.name.append(name)
            dtos.append(dto)

        dtos = annotate_organisms(dtos=dtos, names=scientific_names)

        dfs = [dto.to_df() for dto in dtos]

        self._df = pd.concat(dfs)
        self._df.reset_index(inplace=True, drop=True)

    def integrate(self, *properties):

        if not properties:
            properties = ('ncbi_taxonomy', 'name')

        snapshot = self.node_snapshot()
        latest_identifier = self.node.latest_identifier()
        *_, numeric_identifier = latest_identifier.split('.')
        numeric_identifier = int(numeric_identifier)

        organisms_to_create = []
        create_identifiers = []
        organisms_to_update = []
        update_identifiers = []

        for i, organism in self._df.iterrows():

            to_create = True

            for prop in properties:

                value = organism.get(prop, '')
                snapshot_values = snapshot.loc[:, prop]

                snapshot_mask: pd.Series = snapshot_values == value

                if snapshot_mask.any():

                    organisms_to_update.append(i)
                    protend_id = snapshot[snapshot_mask, self.node.identifying_property]
                    update_identifiers.append(protend_id)

                    to_create = False
                    break

            if to_create:
                organisms_to_create.append(i)
                numeric_identifier += 1
                protend_id = f'{self.node.header}.{self.node.entity}.{numeric_identifier}'
                create_identifiers.append(protend_id)

        organisms_to_create_df = self._df.loc[organisms_to_create, :]
        organisms_to_create_df[self.node.identifying_property] = create_identifiers
        created_nodes = self.node.node_from_df(organisms_to_create_df)

        organisms_to_update_df = self._df.loc[organisms_to_update, :]
        organisms_to_update_df[self.node.identifying_property] = update_identifiers
        updated_nodes = self.node.node_update_from_df(organisms_to_update)