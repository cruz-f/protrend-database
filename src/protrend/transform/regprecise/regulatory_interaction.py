import pandas as pd

from protrend.io import read_json_frame, read_from_stack
from protrend.model import RegulatoryInteraction, Source, Organism, Effector, Regulator, Operon, Gene, TFBS
from protrend.utils.processors import (apply_processors, to_list, to_int_str, regulatory_effect_regprecise,
                                       take_last)
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.regprecise.effector import EffectorTransformer
from protrend.transform.regprecise.operon import OperonTransformer
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.regprecise.source import SourceTransformer
from protrend.utils import SetList


class RegulatoryInteractionTransformer(RegPreciseTransformer,
                                       source='regprecise',
                                       version='0.0.0',
                                       node=RegulatoryInteraction,
                                       order=50,
                                       register=True):
    default_transform_stack = {'effector': 'integrated_effector.json',
                               'regulator': 'integrated_regulator.json',
                               'operon': 'integrated_operon.json'}
    columns = SetList(['regulator', 'regulatory_effect', 'organism_protrend_id', 'url',
                       'regulon_id', 'effector', 'operon', 'genes', 'tfbss',
                       'regulatory_interaction_hash', 'protrend_id'])

    def transform(self) -> pd.DataFrame:
        # merge effector and regulator
        effector = read_from_stack(stack=self.transform_stack, file='effector',
                                   default_columns=EffectorTransformer.columns, reader=read_json_frame)
        effector = self.select_columns(effector, 'protrend_id', 'effector_id')
        effector = effector.rename(columns={'protrend_id': 'effector'})
        effector = apply_processors(effector, effector_id=to_int_str)

        regulator = read_from_stack(stack=self.transform_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = self.select_columns(regulator, 'protrend_id', 'effector', 'operon', 'regulation_mode',
                                        'organism_protrend_id', 'url', 'regulon_id')
        regulator = regulator.rename(columns={'protrend_id': 'regulator',
                                              'effector': 'regulator_effector',
                                              'operon': 'regulator_operon',
                                              'regulation_mode': 'regulatory_effect'})
        regulator = apply_processors(regulator, regulator_effector=to_list, regulator_operon=to_list)
        regulator = regulator.explode(column='regulator_effector')
        regulator = apply_processors(regulator, regulator_effector=to_int_str)

        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = self.select_columns(operon, 'protrend_id', 'operon_id_old', 'genes', 'tfbss')
        operon = operon.rename(columns={'protrend_id': 'operon'})
        operon = apply_processors(operon, operon_id_old=to_list, genes=to_list, tfbss=to_list)
        operon = operon.explode(column='operon_id_old')

        regulator_effector = pd.merge(regulator, effector, how='left',
                                      left_on='regulator_effector', right_on='effector_id')
        regulator_effector = regulator_effector.drop(columns=['effector_id', 'regulator_effector'])
        regulator_effector = regulator_effector.rename(columns={'effector': 'regulator_effector'})

        regulator_effector = apply_processors(regulator_effector, regulator_operon=to_list)
        regulator_effector = regulator_effector.explode(column='regulator_operon')

        regulatory_interaction = pd.merge(regulator_effector, operon,
                                          left_on='regulator_operon', right_on='operon_id_old')
        regulatory_interaction = regulatory_interaction.drop(columns=['operon_id_old', 'regulator_operon'])

        regulatory_interaction = apply_processors(regulatory_interaction,
                                                  regulatory_effect=regulatory_effect_regprecise)

        regulatory_interaction = self.regulatory_interaction_hash(regulatory_interaction)

        self.stack_transformed_nodes(regulatory_interaction)

        return regulatory_interaction


class RegulatoryInteractionToSourceConnector(RegPreciseConnector,
                                             source='regprecise',
                                             version='0.0.0',
                                             from_node=RegulatoryInteraction,
                                             to_node=Source,
                                             register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json',
                             'source': 'integrated_source.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=rin['url'].tolist(),
                      external_identifier=rin['regulon_id'].tolist(),
                      key=['regulon_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class RegulatoryInteractionToOrganismConnector(RegPreciseConnector,
                                               source='regprecise',
                                               version='0.0.0',
                                               from_node=RegulatoryInteraction,
                                               to_node=Organism,
                                               register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        rin = apply_processors(rin, organism_protrend_id=[to_list, take_last])

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['organism_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToEffectorConnector(RegPreciseConnector,
                                               source='regprecise',
                                               version='0.0.0',
                                               from_node=RegulatoryInteraction,
                                               to_node=Effector,
                                               register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, effectors=to_list)
        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['regulator_effector'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToRegulatorConnector(RegPreciseConnector,
                                                source='regprecise',
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


class RegulatoryInteractionToOperonConnector(RegPreciseConnector,
                                             source='regprecise',
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


class RegulatoryInteractionToGeneConnector(RegPreciseConnector,
                                           source='regprecise',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):
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


class RegulatoryInteractionToTFBSConnector(RegPreciseConnector,
                                           source='regprecise',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=TFBS,
                                           register=True):
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
