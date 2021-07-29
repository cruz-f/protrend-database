from typing import Dict, Tuple, List, Callable

import pandas as pd

from protrend.model.model import Source, RegulatoryFamily
from protrend.model.node import protrend_id_decoder, protrend_id_encoder
from protrend.transform.processors import remove_white_space, remove_regprecise_more, remove_multiple_white_space, \
    rstrip, lstrip, remove_pubmed
from protrend.transform.regprecise.settings import RegPreciseTransformSettings
from protrend.transform.transformer import Transformer


class RegulatoryFamilyTransformer(Transformer):
    node = RegulatoryFamily
    integration_properties = ('name', 'rfam')
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
            files = RegPreciseTransformSettings.effector

        super().__init__(source=source, version=version, **files)

    def read(self, *args, **kwargs):
        self.read_json_lines()

    @property
    def tf_family(self):
        df = self.get('tf_family', pd.DataFrame(columns=['regprecise_name']))

        if not df.empty:
            df = df.rename(columns={'name': 'regprecise_name',
                                    'url': 'tf_family_url',
                                    'description': 'tf_family_description',
                                    'pubmed': 'tf_family_pubmed'})

            df = df.drop_duplicates(subset=['regprecise_name'])

            df['regprecise_name'] = df['regprecise_name'].map(remove_white_space)

            df['description'] = df['description'].map(remove_regprecise_more)
            df['description'] = df['description'].map(remove_pubmed)

            df['description'] = df['description'].map(remove_multiple_white_space)
            df['description'] = df['description'].map(rstrip)
            df['description'] = df['description'].map(lstrip)

        return df

    @property
    def tf(self):
        df = self.get('tf', pd.DataFrame(columns=['regprecise_name']))

        if not df.empty:
            df = df.rename(columns={'name': 'regprecise_name',
                                    'url': 'tf_url',
                                    'description': 'tf_description',
                                    'pubmed': 'tf_pubmed'})

            df = df.drop_duplicates(subset=['regprecise_name'])

            df['regprecise_name'] = df['regprecise_name'].map(remove_white_space)

            df['description'] = df['description'].map(remove_regprecise_more)
            df['description'] = df['description'].map(remove_pubmed)

            df['description'] = df['description'].map(remove_multiple_white_space)
            df['description'] = df['description'].map(rstrip)
            df['description'] = df['description'].map(lstrip)

        return df

    @property
    def rna(self):
        df = self.get('rna', pd.DataFrame(columns=['regprecise_name']))

        if not df.empty:
            df = df.rename(columns={'name': 'regprecise_name',
                                    'url': 'rna_url',
                                    'description': 'rna_description',
                                    'pubmed': 'rna_pubmed'})

            df = df.drop_duplicates(subset=['regprecise_name'])

            df['regprecise_name'] = df['regprecise_name'].map(remove_white_space)

            df['description'] = df['description'].map(remove_regprecise_more)
            df['description'] = df['description'].map(remove_pubmed)

            df['description'] = df['description'].map(remove_multiple_white_space)
            df['description'] = df['description'].map(rstrip)
            df['description'] = df['description'].map(lstrip)

        return df

    def transform(self) -> pd.DataFrame:

        df = pd.merge(self.tf_family, self.tf, on='regprecise_name')
        df = pd.merge(df, self.rna, on='regprecise_name')

        # TODO: missing merge of the multiple columns into one

        return df

    @staticmethod
    def find_regulatory_family(regulatory_family: pd.Series, snapshot_properties: Dict[str, pd.Series]) -> pd.Series:

        for prop, snapshot_values in snapshot_properties.items():

            identifier = regulatory_family.get(prop, None)

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

        # TODO: missing

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
        for i, family in df.iterrows():

            family_mask = self.find_regulatory_family(family, snapshot_properties)

            # update
            if family_mask.any():
                protend_id = snapshot.loc[family_mask, self.node.identifying_property].iloc[0]
                to_update.append((i, protend_id))

            else:
                integer += 1
                protend_id = protrend_id_encoder(self.node.header, self.node.entity, integer)
                to_create.append((i, protend_id))

        to_create = self.integrate_nodes(df=df, index=to_create, node_factory=self.node.node_from_df)
        to_update = self.integrate_nodes(df=df, index=to_update, node_factory=self.node.node_update_from_df)

        df = pd.concat([to_create, to_update])
        self.stack_csv('regulatory_family', df)
        return df

    def load_relationships(self, df: pd.DataFrame) -> Dict[str, pd.DataFrame]:

        # TODO: missing

        n_rows, _ = df.shape

        from_ = ['regulatory_family'] * n_rows
        to = ['source'] * n_rows
        from_property = ['protrend_id'] * n_rows
        to_property = ['name'] * n_rows

        source_df = pd.DataFrame([from_, to, from_property, to_property],
                                 columns=['from', 'to', 'from_property', 'to_property'])

        source_df['name'] = ['effector_id'] * n_rows
        source_df['url'] = df.loc[:, 'url']
        source_df['external_identifier'] = df.loc[:, 'effector_id']

        self.stack_csv('regulatory_family_source', source_df)

        return {(RegulatoryFamily.node_name(), Source.node_name()): df}
