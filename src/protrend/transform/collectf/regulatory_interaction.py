import pandas as pd

from protrend.io.json import read_json_frame, read_json_lines
from protrend.io.utils import read_from_stack
from protrend.model.model import RegulatoryInteraction, Regulator, Operon, Gene, TFBS
from protrend.transform.collectf.base import CollectfTransformer, CollectfConnector
from protrend.transform.collectf.operon import OperonTransformer
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.transform.collectf.tfbs import TFBSTransformer
from protrend.transform.processors import (apply_processors, to_list, regulatory_effect_collectf, to_set, flatten_set,
                                           take_first, to_list_nan)


class RegulatoryInteractionTransformer(CollectfTransformer):
    default_node = RegulatoryInteraction
    default_node_factors = ()
    default_transform_stack = {'regulator': 'integrated_regulator.json',
                               'operon': 'integrated_operon.json',
                               'gene': 'integrated_gene.json',
                               'tfbs': 'TFBS.json'}
    default_order = 50
    columns = {'protrend_id',
               'regulator', 'operon', 'genes', 'tfbss', 'regulatory_effect', 'organism_protrend_id',
               'regulon'}

    def transform(self) -> pd.DataFrame:
        tfbs = read_from_stack(stack=self.transform_stack, file='tfbs',
                               default_columns=TFBSTransformer.read_columns, reader=read_json_lines)
        tfbs = self.select_columns(tfbs, 'regulon', 'operon', 'mode')
        tfbs = apply_processors(tfbs, regulon=to_list_nan, operon=to_list_nan)
        tfbs = tfbs.explode('regulon')
        tfbs = tfbs.explode('operon')
        tfbs = tfbs.dropna(subset=['regulon'])
        tfbs = tfbs.dropna(subset=['operon'])
        tfbs = tfbs.drop_duplicates(subset=['regulon', 'operon'])

        regulator = read_from_stack(stack=self.transform_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = self.select_columns(regulator, 'protrend_id', 'uniprot_accession', 'organism_protrend_id')
        regulator = regulator.rename(columns={'protrend_id': 'regulator', 'uniprot_accession': 'regulon'})

        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = self.select_columns(operon, 'protrend_id', 'genes', 'tfbss', 'operon_id_old')
        operon = operon.rename(columns={'operon_id_old': 'operon', 'protrend_id': 'operon_protrend_id'})
        operon = apply_processors(operon, genes=to_list_nan, tfbss=to_list_nan, operon=to_list_nan)
        operon = operon.explode(column='operon')

        df = pd.merge(tfbs, regulator, on='regulon')
        df = df.dropna(subset=['regulator'])

        df = pd.merge(df, operon, on='operon')
        df = df.dropna(subset=['operon_protrend_id'])

        df = df.rename(columns={'operon': 'operon_id', 'operon_protrend_id': 'operon', 'mode': 'regulatory_effect'})
        df = self.drop_duplicates(df=df, subset=['regulator', 'operon'], perfect_match=True, preserve_nan=True)
        df = df.dropna(subset=['regulator'])
        df = df.dropna(subset=['operon'])

        aggregation = {'genes': flatten_set, 'tfbss': flatten_set,
                       'organism_protrend_id': to_set, 'regulator': to_set}
        df = self.group_by(df=df, column='operon', aggregation=aggregation, default=take_first)
        mask = df['organism_protrend_id'].map(len) == 1
        df = df[mask]
        df = apply_processors(df, organism_protrend_id=[to_list, take_first], regulator=to_list)
        df = df.explode(column='regulator')

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


class RegulatoryInteractionToRegulatorConnector(CollectfConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Regulator
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToOperonConnector(CollectfConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Operon
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToGeneConnector(CollectfConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Gene
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, genes=to_list)
        rin = rin.explode(column='genes')
        rin = rin.dropna(subset=['genes'])

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['genes'].tolist()

        kwargs = dict(operon=rin['operon'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatoryInteractionToTFBSConnector(CollectfConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = TFBS
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, tfbss=to_list)
        rin = rin.explode(column='tfbss')
        rin = rin.dropna(subset=['tfbss'])

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['tfbss'].tolist()

        kwargs = dict(operon=rin['operon'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatorToOperonConnector(CollectfConnector):
    default_from_node = Regulator
    default_to_node = Operon
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatorToGeneConnector(CollectfConnector):
    default_from_node = Regulator
    default_to_node = Gene
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, genes=to_list)
        rin = rin.explode(column='genes')

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['genes'].tolist()
        kwargs = dict(operon=rin['operon'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatorToTFBSConnector(CollectfConnector):
    default_from_node = Regulator
    default_to_node = TFBS
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, tfbss=to_list)
        rin = rin.explode(column='tfbss')

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['tfbss'].tolist()
        kwargs = dict(operon=rin['operon'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)