import pandas as pd

from protrend.io import read_json_frame, read_from_stack
from protrend.model import RegulatoryInteraction, Regulator, Operon, Gene, TFBS
from protrend.transform.dbtbs.base import DBTBSTransformer, DBTBSConnector
from protrend.transform.dbtbs.operon import OperonTransformer
from protrend.transform.dbtbs.regulator import RegulatorTransformer
from protrend.transform.dbtbs.tfbs import TFBSTransformer
from protrend.utils.processors import (apply_processors, take_first, regulatory_effect_dbtbs, to_list)
from protrend.utils import SetList


class RegulatoryInteractionTransformer(DBTBSTransformer,
                                       source='dbtbs',
                                       version='0.0.3',
                                       node=RegulatoryInteraction,
                                       order=80,
                                       register=True):
    default_transform_stack = {'regulator': 'integrated_regulator.json',
                               'operon': 'integrated_operon.json',
                               'tfbs': 'integrated_tfbs.json'}
    columns = SetList(['identifier', 'url', 'pubmed', 'operon_name',
                       'regulator', 'operon', 'genes', 'tfbss', 'regulator_effector', 'regulatory_effect'])

    def _transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = self.select_columns(tfbs, 'identifier', 'url', 'regulation', 'pubmed', 'tf', 'operon')
        tfbs = tfbs.explode(column='tf')
        tfbs = tfbs.explode(column='operon')
        tfbs = self.drop_duplicates(df=tfbs, subset=['tf', 'operon'], perfect_match=True, preserve_nan=True)
        tfbs = tfbs.dropna(subset=['tf', 'operon'])
        return tfbs

    def _transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = self.select_columns(regulator, 'protrend_id', 'name_dbtbs')
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id',
                                              'name_dbtbs': 'regulator_name_dbtbs'})
        return regulator

    def _transform_operon(self, operon: pd.DataFrame) -> pd.DataFrame:
        operon = self.select_columns(operon, 'protrend_id', 'name', 'genes', 'tfbss')
        operon = operon.rename(columns={'protrend_id': 'operon_protrend_id', 'name': 'operon_name'})
        return operon

    def transform(self) -> pd.DataFrame:
        tfbs = read_from_stack(stack=self.transform_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = self._transform_tfbs(tfbs)

        regulator = read_from_stack(stack=self.transform_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = self._transform_regulator(regulator)

        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = self._transform_operon(operon)

        regulator_tfbs = pd.merge(regulator, tfbs, left_on='regulator_name_dbtbs', right_on='tf')

        regulatory_interaction = pd.merge(regulator_tfbs, operon, left_on='operon', right_on='operon_name')

        regulatory_interaction = regulatory_interaction.drop(columns=['regulator_name_dbtbs', 'operon', 'tf'])
        regulatory_interaction = regulatory_interaction.rename(columns={'regulator_protrend_id': 'regulator',
                                                                        'operon_protrend_id': 'operon',
                                                                        'regulation': 'regulatory_effect'})

        regulatory_interaction = apply_processors(regulatory_interaction,
                                                  regulatory_effect=[take_first, regulatory_effect_dbtbs])
        regulatory_interaction['regulator_effector'] = None

        regulatory_interaction = self.regulatory_interaction_hash(regulatory_interaction)

        self._stack_transformed_nodes(regulatory_interaction)

        return regulatory_interaction


class RegulatoryInteractionToRegulatorConnector(DBTBSConnector,
                                                source='dbtbs',
                                                version='0.0.3',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToOperonConnector(DBTBSConnector,
                                             source='dbtbs',
                                             version='0.0.3',
                                             from_node=RegulatoryInteraction,
                                             to_node=Operon,
                                             register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToGeneConnector(DBTBSConnector,
                                           source='dbtbs',
                                           version='0.0.3',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):
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


class RegulatoryInteractionToTFBSConnector(DBTBSConnector,
                                           source='dbtbs',
                                           version='0.0.3',
                                           from_node=RegulatoryInteraction,
                                           to_node=TFBS,
                                           register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
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


class RegulatorToOperonConnector(DBTBSConnector,
                                 source='dbtbs',
                                 version='0.0.3',
                                 from_node=Regulator,
                                 to_node=Operon,
                                 register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatorToGeneConnector(DBTBSConnector,
                               source='dbtbs',
                               version='0.0.3',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):
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


class RegulatorToTFBSConnector(DBTBSConnector,
                               source='dbtbs',
                               version='0.0.3',
                               from_node=Regulator,
                               to_node=TFBS,
                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
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
