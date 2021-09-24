import pandas as pd

from protrend.io.json import read_json_frame
from protrend.io.utils import read_from_stack
from protrend.model.model import RegulatoryInteraction, Source, Organism, Effector, Regulator, Operon, Gene, TFBS
from protrend.transform.processors import (apply_processors, to_list, to_int_str, to_set_list, flatten_set_list,
                                           regulatory_effect, take_last, regulatory_interaction_hash)
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.regprecise.effector import EffectorTransformer
from protrend.transform.regprecise.operon import OperonTransformer
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.regprecise.source import SourceTransformer
from protrend.utils import SetList


class RegulatoryInteractionTransformer(RegPreciseTransformer):
    default_node = RegulatoryInteraction
    default_node_factors = SetList(['regulatory_interaction_hash'])
    default_transform_stack = {'effector': 'integrated_effector.json',
                               'regulator': 'integrated_regulator.json',
                               'operon': 'integrated_operon.json'}
    default_order = 50
    columns = SetList(['regulator', 'regulatory_effect', 'organism_protrend_id', 'url',
                       'regulon_id', 'effectors', 'operon', 'genes', 'tfbss',
                       'regulatory_interaction_hash', 'protrend_id'])

    def transform(self) -> pd.DataFrame:
        # merge effector and regulator
        effector = read_from_stack(stack=self.transform_stack, file='effector',
                                   default_columns=EffectorTransformer.columns, reader=read_json_frame)
        effector = self.select_columns(effector, 'protrend_id', 'effector_id')
        effector = effector.rename(columns={'protrend_id': 'effectors'})
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
        df = pd.merge(regulator, effector, how='left', left_on='regulator_effector', right_on='effector_id')

        aggregation = {'regulator_operon': flatten_set_list,
                       'organism_protrend_id': take_last, 'url': take_last, 'regulon_id': take_last}
        df = self.group_by(df=df, column='regulator', aggregation=aggregation, default=to_set_list)
        df = df.drop(columns=['effector_id', 'regulator_effector'])

        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = self.select_columns(operon, 'protrend_id', 'operon_id_old', 'genes', 'tfbss')
        operon = operon.rename(columns={'protrend_id': 'operon'})
        operon = apply_processors(operon, operon_id_old=to_list, genes=to_list, tfbss=to_list)

        df = apply_processors(df, regulator_operon=to_list)
        df = df.explode(column='regulator_operon')
        operon = operon.explode(column='operon_id_old')
        df = pd.merge(df, operon, left_on='regulator_operon', right_on='operon_id_old')
        df = df.drop(columns=['operon_id_old', 'regulator_operon'])
        df = self.drop_duplicates(df=df, subset=['regulator', 'operon'], perfect_match=True, preserve_nan=True)
        df = df.dropna(subset=['regulator', 'operon'])

        # filter by regulator + operon
        df2 = apply_processors(df, regulator=to_list, operon=to_list)

        df['regulatory_interaction_hash'] = df2['regulator'] + df2['operon']
        df = apply_processors(df, regulatory_interaction_hash=regulatory_interaction_hash)
        df = df.dropna(subset=['regulatory_interaction_hash'])

        df = apply_processors(df, regulatory_effect=regulatory_effect)

        self._stack_transformed_nodes(df)

        return df


class RegulatoryInteractionToSourceConnector(RegPreciseConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Source
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


class RegulatoryInteractionToOrganismConnector(RegPreciseConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Organism
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


class RegulatoryInteractionToEffectorConnector(RegPreciseConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Effector
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        rin = apply_processors(rin, effectors=to_list)
        rin = rin.explode(column='effectors')
        rin = rin.dropna(subset=['effectors'])
        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['effectors'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToRegulatorConnector(RegPreciseConnector):
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


class RegulatoryInteractionToOperonConnector(RegPreciseConnector):
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


class RegulatoryInteractionToGeneConnector(RegPreciseConnector):
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


class RegulatoryInteractionToTFBSConnector(RegPreciseConnector):
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
