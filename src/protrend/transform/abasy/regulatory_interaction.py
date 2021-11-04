import pandas as pd

from protrend.io.json import read_json_frame
from protrend.io.utils import read_from_stack
from protrend.model.model import RegulatoryInteraction, Regulator, Operon, Gene
from protrend.transform.abasy.base import AbasyTransformer, AbasyConnector
from protrend.transform.abasy.operon import OperonTransformer
from protrend.transform.abasy.regulator import RegulatorTransformer
from protrend.utils.processors import (apply_processors, to_list, regulatory_effect_coryneregnet, rstrip, lstrip)
from protrend.utils import SetList


class RegulatoryInteractionTransformer(AbasyTransformer,
                                       source='abasy',
                                       version='0.0.0',
                                       node=RegulatoryInteraction,
                                       order=80,
                                       register=True):
    default_transform_stack = {'regulator': 'integrated_regulator.json',
                               'operon': 'integrated_operon.json'}
    columns = SetList(['regulator_effector', 'regulator', 'operon', 'genes', 'tfbss', 'regulatory_effect',
                       'regulatory_interaction_hash', 'protrend_id',
                       'id', 'source', 'target', 'Effect', 'Evidence', 'taxonomy', 'source_target_taxonomy'])

    def _transform_networks(self, networks: pd.DataFrame) -> pd.DataFrame:
        networks = networks.dropna(subset=['source', 'target', 'taxonomy', 'Effect'])
        net_src_t_tax = networks['source'] + networks['target'] + networks['Effect'] + networks['taxonomy']
        networks.loc[:, 'source_target_taxonomy'] = net_src_t_tax

        networks = self.drop_duplicates(df=networks, subset=['source_target_taxonomy'], perfect_match=True)
        networks = networks.dropna(subset=['source_target_taxonomy'])

        networks = apply_processors(networks,
                                    source_target_taxonomy=[rstrip, lstrip],
                                    source=[rstrip, lstrip],
                                    target=[rstrip, lstrip])

        return networks

    def _transform_regulator(self) -> pd.DataFrame:
        regulator = read_from_stack(stack=self.transform_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = self.select_columns(regulator, 'protrend_id', 'source')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def _transform_operon(self) -> pd.DataFrame:
        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = self.select_columns(operon, 'protrend_id', 'genes', 'tfbss', 'target')
        operon = operon.rename(columns={'protrend_id': 'operon'})
        operon = apply_processors(operon, genes=to_list, tfbss=to_list)
        return operon

    def transform(self) -> pd.DataFrame:
        networks = self._build_networks()
        network = self._transform_networks(networks)

        regulator = self._transform_regulator()
        regulator_network = pd.merge(network, regulator, on='source')

        operon = self._transform_operon()

        regulatory_interaction = pd.merge(regulator_network, operon, on='target')
        regulatory_interaction.loc[:, 'regulatory_effect'] = regulatory_interaction['Effect'].tolist()

        regulatory_interaction = apply_processors(regulatory_interaction,
                                                  regulatory_effect=regulatory_effect_coryneregnet)

        regulatory_interaction['regulator_effector'] = None

        regulatory_interaction = self.regulatory_interaction_hash(regulatory_interaction)

        self._stack_transformed_nodes(regulatory_interaction)

        return regulatory_interaction


class RegulatoryInteractionToRegulatorConnector(AbasyConnector,
                                                source='abasy',
                                                version='0.0.0',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=False):

    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToOperonConnector(AbasyConnector,
                                             source='abasy',
                                             version='0.0.0',
                                             from_node=RegulatoryInteraction,
                                             to_node=Operon,
                                             register=False):

    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToGeneConnector(AbasyConnector,
                                           source='abasy',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=False):

    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
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


class RegulatorToOperonConnector(AbasyConnector,
                                 source='abasy',
                                 version='0.0.0',
                                 from_node=Regulator,
                                 to_node=Operon,
                                 register=False):

    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatorToGeneConnector(AbasyConnector,
                               source='abasy',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=Gene,
                               register=False):

    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
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
