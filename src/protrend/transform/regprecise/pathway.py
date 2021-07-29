from typing import Dict, Tuple, List, Callable

import pandas as pd

from protrend.model.model import Source, Pathway
from protrend.model.node import protrend_id_decoder, protrend_id_encoder
from protrend.transform.annotation.pathway import annotate_pathways
from protrend.transform.dto import PathwayDTO
from protrend.transform.processors import rstrip, lstrip
from protrend.transform.regprecise.settings import RegPreciseTransformSettings
from protrend.transform.transformer import Transformer


class PathwayTransformer(Transformer):
    node = Pathway
    integration_properties = ('name', )
    source_name = 'regprecise'

    def __init__(self,
                 source: str = None,
                 version: str = None,
                 **files: Dict[str, str]):

        if not source:
            source = RegPreciseTransformSettings.source

        if not version:
            version = RegPreciseTransformSettings.version

        if not files:
            files = RegPreciseTransformSettings.pathway

        super().__init__(source=source, version=version, **files)

    def read(self, *args, **kwargs):
        self.read_json_lines()

    @property
    def pathway(self):
        df = self.get('pathway')

        if df is not None:
            df = df.rename(columns={'name': 'regprecise_name'})
            df = df.drop_duplicates(subset=['regprecise_name'])
            df['regprecise_name'] = df['regprecise_name'].map(rstrip)
            df['regprecise_name'] = df['regprecise_name'].map(lstrip)

        return df

    def transform(self) -> pd.DataFrame:

        if self.pathway is None:
            return pd.DataFrame()

        names = list(self.pathway.loc[:, 'regprecise_name'])

        dtos = [PathwayDTO(input_value=name) for name in names]
        annotate_pathways(dtos=dtos, names=names)

        pathways = pd.DataFrame([dto.to_dict() for dto in dtos])

        df = pd.merge(self.pathway, pathways, left_on='regprecise_name', right_on='input_value')

        return df

    @staticmethod
    def find_pathway(pathway: pd.Series, snapshot_properties: Dict[str, pd.Series]) -> pd.Series:

        for prop, snapshot_values in snapshot_properties.items():

            identifier = pathway.get(prop, None)

            if identifier is None:
                continue

            snapshot_mask: pd.Series = snapshot_values == identifier

            if snapshot_mask.any():
                return snapshot_mask

        return pd.Series([False])

    def integrate_nodes(self,
                        df: pd.DataFrame,
                        index: List[Tuple[int, str]],
                        node_factory: Callable):

        to_idx, to_ids = list(zip(*index))
        to_df = df.loc[to_ids, :]
        to_df[self.node.identifying_property] = to_idx

        node_factory(nodes=to_df, save=True)
        return to_df

    def load_nodes(self, df, *properties) -> pd.DataFrame:

        snapshot = self.node_snapshot()

        if not properties:
            properties = self.integration_properties

        snapshot_properties = {prop: snapshot.loc[:, prop] for prop in properties}

        last_node = self.node.last_node()
        if last_node is None:
            integer = 0

        else:
            integer = protrend_id_decoder(last_node.protrend_id)

        to_update: List[Tuple[int, str]] = []
        to_create: List[Tuple[int, str]] = []
        for i, pathway in df.iterrows():

            pathway_mask = self.find_pathway(pathway, snapshot_properties)

            # update
            if pathway_mask.any():
                protend_id = snapshot.loc[pathway_mask, self.node.identifying_property].iloc[0]
                to_update.append((i, protend_id))

            else:
                integer += 1
                protend_id = protrend_id_encoder(self.node.header, self.node.entity, integer)
                to_create.append((i, protend_id))

        to_create = self.integrate_nodes(df=df, index=to_create, node_factory=self.node.node_from_df)
        to_update = self.integrate_nodes(df=df, index=to_update, node_factory=self.node.node_update_from_df)

        df = pd.concat([to_create, to_update])
        self.stack_csv('pathway', df)
        return df

    def load_relationships(self, df: pd.DataFrame) -> Dict[str, pd.DataFrame]:

        n_rows, _ = df.shape

        from_ = ['pathway'] * n_rows
        to = ['source'] * n_rows
        from_property = ['protrend_id'] * n_rows
        to_property = ['name'] * n_rows

        source_df = pd.DataFrame([from_, to, from_property, to_property],
                                 columns=['from', 'to', 'from_property', 'to_property'])

        source_df['name'] = ['pathway_id'] * n_rows
        source_df['url'] = df.loc[:, 'url']
        source_df['external_identifier'] = df.loc[:, 'pathway_id']

        self.stack_csv('pathway_source', source_df)

        return {(Pathway.node_name(), Source.node_name()): df}
