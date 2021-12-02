import pandas as pd

from protrend.io import read_json_frame, read_from_stack
from protrend.model import RegulatoryInteraction, Effector, Regulator, Operon, Gene
from protrend.transform.literature.base import LiteratureTransformer, LiteratureConnector
from protrend.transform.literature.effector import EffectorTransformer
from protrend.transform.literature.operon import OperonTransformer
from protrend.transform.literature.regulator import RegulatorTransformer
from protrend.utils.processors import (apply_processors, to_list, regulatory_effect_literature, to_set_list)
from protrend.utils import SetList


class RegulatoryInteractionTransformer(LiteratureTransformer,
                                       source='literature',
                                       version='0.0.0',
                                       node=RegulatoryInteraction,
                                       order=80,
                                       register=True):
    default_transform_stack = {'effector': 'integrated_effector.json',
                               'regulator': 'integrated_regulator.json',
                               'operon': 'integrated_operon.json'}
    columns = SetList(['regulator_effector', 'regulator', 'operon', 'genes', 'tfbss', 'regulatory_effect',
                       'regulatory_interaction_hash', 'protrend_id',
                       'regulator_locus_tag', 'operon', 'genes_locus_tag',
                       'regulatory_effect', 'evidence', 'effector', 'mechanism',
                       'publication', 'taxonomy', 'source', 'network_id'])

    def _transform_effector(self) -> pd.DataFrame:
        effector = read_from_stack(stack=self.transform_stack, file='effector',
                                   default_columns=EffectorTransformer.columns, reader=read_json_frame)
        effector = self.select_columns(effector, 'protrend_id', 'network_id')
        effector = effector.rename(columns={'protrend_id': 'regulator_effector'})
        return effector

    def _transform_regulator(self) -> pd.DataFrame:
        regulator = read_from_stack(stack=self.transform_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = self.select_columns(regulator, 'protrend_id', 'regulator_locus_tag')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def _transform_operon(self) -> pd.DataFrame:
        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = self.select_columns(operon, 'protrend_id', 'genes', 'operon_id')
        operon = operon.rename(columns={'protrend_id': 'operon'})
        operon = apply_processors(operon, genes=to_list)
        return operon

    def transform(self) -> pd.DataFrame:
        # 'regulator_locus_tag', 'operon', 'genes_locus_tag',
        # 'regulatory_effect', 'evidence', 'effector', 'mechanism',
        # 'publication', 'taxonomy', 'source', 'network_id'
        network = self._build_network()
        network['operon_id'] = network['operon'] + network['taxonomy']
        network = network.drop(columns=['operon'])

        # 'regulator_effector', 'network_id'
        effector = self._transform_effector()

        # 'regulator', 'network_id'
        regulator = self._transform_regulator()

        # 'operon', 'genes', 'network_id'
        operon = self._transform_operon()

        # 'regulator_locus_tag', 'operon', 'genes_locus_tag',
        # 'regulatory_effect', 'evidence', 'effector', 'mechanism',
        # 'publication', 'taxonomy', 'source', 'network_id', 'regulator_effector'
        network_effector = pd.merge(network, effector, how='left', on='network_id')

        # 'regulator_locus_tag', 'operon', 'genes_locus_tag',
        # 'regulatory_effect', 'evidence', 'effector', 'mechanism',
        # 'publication', 'taxonomy', 'source', 'network_id', 'regulator_effector', 'regulator'
        network_effector_regulator = pd.merge(network_effector, regulator, on='regulator_locus_tag')

        # 'regulator_locus_tag', 'operon', 'genes_locus_tag',
        # 'regulatory_effect', 'evidence', 'effector', 'mechanism',
        # 'publication', 'taxonomy', 'source', 'network_id', 'regulator_effector', 'regulator', 'operon', 'genes'
        regulatory_interaction = pd.merge(network_effector_regulator, operon, on='operon_id')

        regulatory_interaction = apply_processors(regulatory_interaction,
                                                  regulatory_effect=regulatory_effect_literature)

        regulatory_interaction = self.regulatory_interaction_hash(regulatory_interaction)

        self.stack_transformed_nodes(regulatory_interaction)

        return regulatory_interaction


class RegulatoryInteractionToEffectorConnector(LiteratureConnector,
                                               source='literature',
                                               version='0.0.0',
                                               from_node=RegulatoryInteraction,
                                               to_node=Effector,
                                               register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, regulator_effector=to_list)
        rin = rin.explode(column='regulator_effector')
        rin = rin.dropna(subset=['regulator_effector'])
        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['regulator_effector'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToRegulatorConnector(LiteratureConnector,
                                                source='literature',
                                                version='0.0.0',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToOperonConnector(LiteratureConnector,
                                             source='literature',
                                             version='0.0.0',
                                             from_node=RegulatoryInteraction,
                                             to_node=Operon,
                                             register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToGeneConnector(LiteratureConnector,
                                           source='literature',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, genes=to_set_list)
        rin = rin.explode(column='genes')
        rin = rin.dropna(subset=['genes'])

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['genes'].tolist()

        kwargs = dict(operon=rin['operon'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatorToEffectorConnector(LiteratureConnector,
                                   source='literature',
                                   version='0.0.0',
                                   from_node=Regulator,
                                   to_node=Effector,
                                   register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, regulator_effector=to_set_list)
        rin = rin.explode(column='regulator_effector')
        rin = rin.dropna(subset=['regulator_effector'])
        rin = rin.drop_duplicates(subset=['regulator', 'regulator_effector'])

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['regulator_effector'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatorToOperonConnector(LiteratureConnector,
                                 source='literature',
                                 version='0.0.0',
                                 from_node=Regulator,
                                 to_node=Operon,
                                 register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatorToGeneConnector(LiteratureConnector,
                               source='literature',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):
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
