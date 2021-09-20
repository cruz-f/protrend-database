import pandas as pd

from protrend.io.json import read_json_frame
from protrend.io.utils import read_from_stack
from protrend.model.model import RegulatoryInteraction
from protrend.transform.collectf.base import CollectfTransformer
from protrend.transform.collectf.operon import OperonTransformer
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.transform.collectf.tfbs import TFBSTransformer
from protrend.transform.processors import (apply_processors, to_list, regulatory_effect_collectf)


class RegulatoryInteractionTransformer(CollectfTransformer):
    default_node = RegulatoryInteraction
    default_node_factors = ()
    default_transform_stack = {'regulator': 'integrated_regulator.json',
                               'operon': 'integrated_operon.json',
                               'gene': 'integrated_gene.json',
                               'tfbs': 'integrated_tfbs.json'}
    default_order = 50
    columns = {'protrend_id',
               'regulator', 'operon', 'genes', 'tfbss', 'regulatory_effect',
               'regulon'}

    def transform(self) -> pd.DataFrame:
        regulator = read_from_stack(stack=self.transform_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = self.select_columns(regulator, 'protrend_id', 'organism_protrend_id', 'uniprot_accession')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})

        tfbs = read_from_stack(stack=self.transform_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = self.select_columns(tfbs, 'regulon', 'operon', 'mode')
        tfbs = apply_processors(tfbs, regulon=to_list)
        tfbs = tfbs.explode(column='regulon')

        df = pd.merge(regulator, tfbs, left_on='uniprot_accession', right_on='regulon')
        df = apply_processors(df, operon=to_list)
        df = df.explode(column='operon')

        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = self.select_columns(operon, 'protrend_id', 'operon_id_old', 'genes', 'tfbss')
        operon = apply_processors(operon, operon_id_old=to_list, genes=to_list, tfbss=to_list)
        operon = operon.explode(column='operon_id_old')

        df = pd.merge(df, operon, left_on='operon', right_on='operon_id_old')

        df = df.drop(columns=['operon_id_old', 'uniprot_accession', 'operon'])
        df = df.rename(columns={'protrend_id': 'operon', 'mode': 'regulatory_effect'})
        df = self.drop_duplicates(df=df, subset=['regulator', 'operon'], perfect_match=True, preserve_nan=True)
        df = df.dropna(subset=['regulator', 'operon'])

        df = apply_processors(df, regulatory_effect=regulatory_effect_collectf)

        self._stack_transformed_nodes(df)

        return df

    def _update_nodes(self, df: pd.DataFrame, mask: pd.Series, snapshot: pd.DataFrame) -> pd.DataFrame:
        # nodes to be updated
        nodes = df[mask]

        if nodes.empty:
            nodes['protrend_id'] = None
            nodes['load'] = None
            nodes['what'] = None
            return nodes

        # find/set protrend identifiers for update nodes
        ids_mask = self.find_snapshot(nodes=nodes, snapshot=snapshot, node_factors=('regulatory_interaction_id',))
        nodes.loc[:, 'protrend_id'] = snapshot.loc[ids_mask, 'protrend_id']
        nodes.loc[:, 'load'] = 'update'
        nodes.loc[:, 'what'] = 'nodes'

        return nodes

    def integrate(self, df: pd.DataFrame) -> pd.DataFrame:
        df['regulatory_interaction_id'] = df['regulator'] + df['operon']

        # ensure uniqueness
        df = self.drop_duplicates(df=df, subset=['regulatory_interaction_id'], perfect_match=True, preserve_nan=True)

        # take a db snapshot for the current node
        snapshot = self.node_view()
        snapshot['regulatory_interaction_id'] = snapshot['regulator'] + snapshot['operon']

        # find matching nodes according to several node factors/properties
        mask = self.find_nodes(nodes=df, snapshot=snapshot, node_factors=('regulatory_interaction_id',))

        # nodes to be updated
        update_nodes = self._update_nodes(df=df, mask=mask, snapshot=snapshot)

        # nodes to be created
        create_nodes = self._create_nodes(df=df, mask=mask)

        # concat both dataframes
        df = pd.concat([create_nodes, update_nodes], axis=0)
        df = df.drop(columns=['regulatory_interaction_id'])

        self._stack_integrated_nodes(df)
        self._stack_nodes(df)

        return df
