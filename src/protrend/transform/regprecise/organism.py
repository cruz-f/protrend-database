from typing import Dict

import pandas as pd

from protrend.model.model import Organism
from protrend.model.node import protrend_id_decoder, protrend_id_encoder
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
        genome: pd.DataFrame = self.get('genome', pd.DataFrame(columns=['name']))

        names = genome.loc[:, 'name']

        dtos = []
        for name in names:
            dto = OrganismDTO()
            dto.name.append(name)
            dtos.append(dto)

        annotate_organisms(dtos=dtos, names=names)

        annotated_df = pd.DataFrame([dto.to_dict() for dto in dtos])

        df = pd.merge(annotated_df, genome, on='name')
        self._df = df

    def integrate(self, *properties):

        if not properties:
            properties = ('ncbi_taxonomy', 'name')

        snapshot = self.node_snapshot()
        latest_identifier = self.node.latest_identifier()
        integer = protrend_id_decoder(latest_identifier)

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
                    protend_id = snapshot.loc[snapshot_mask, self.node.identifying_property].iloc[0]
                    update_identifiers.append(protend_id)

                    to_create = False
                    break

            if to_create:
                organisms_to_create.append(i)
                integer += 1
                protend_id = protrend_id_encoder(self.node.header, self.node.entity, integer)
                create_identifiers.append(protend_id)

        organisms_to_create_df = self._df.loc[organisms_to_create, :]
        organisms_to_create_df[self.node.identifying_property] = create_identifiers
        self.node.node_from_df(organisms_to_create_df, save=True)

        organisms_to_update_df = self._df.loc[organisms_to_update, :]
        organisms_to_update_df[self.node.identifying_property] = update_identifiers
        self.node.node_update_from_df(organisms_to_update_df, save=True)

        df = pd.concat([organisms_to_create_df, organisms_to_update_df])

        self.stack_csv('organism', df)
