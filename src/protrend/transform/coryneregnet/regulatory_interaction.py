import pandas as pd

from protrend.io import read_json_frame, read_from_stack
from protrend.model import RegulatoryInteraction, Regulator, Operon, Gene, TFBS
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer, CoryneRegNetConnector
from protrend.transform.coryneregnet.operon import OperonTransformer
from protrend.transform.coryneregnet.regulator import RegulatorTransformer
from protrend.utils.processors import (apply_processors, to_list, regulatory_effect_coryneregnet)
from protrend.utils import SetList


class RegulatoryInteractionTransformer(CoryneRegNetTransformer,
                                       source='coryneregnet',
                                       version='0.0.0',
                                       node=RegulatoryInteraction,
                                       order=70,
                                       register=True):
    default_transform_stack = {'regulator': 'integrated_regulator.json',
                               'operon': 'integrated_operon.json'}
    columns = SetList(['regulator_effector', 'regulator', 'operon', 'genes', 'tfbss', 'regulatory_effect',
                       'regulatory_interaction_hash', 'protrend_id',
                       'Operon', 'Orientation', 'Genes', 'TF_locusTag', 'TG_locusTag',
                       'Role', 'Evidence', 'PMID', 'Source', 'taxonomy'])

    def _transform_regulator(self) -> pd.DataFrame:
        regulator = read_from_stack(stack=self.transform_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = self.select_columns(regulator, 'protrend_id', 'TF_locusTag')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def _transform_operon(self) -> pd.DataFrame:
        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = self.select_columns(operon, 'protrend_id', 'genes', 'tfbss', 'Operon', 'Orientation', 'Genes')
        operon = operon.rename(columns={'protrend_id': 'operon'})
        operon = apply_processors(operon, genes=to_list, tfbss=to_list, Genes=to_list)
        return operon

    def transform(self) -> pd.DataFrame:
        # 'TF_locusTag', 'TG_locusTag', 'Role', 'Evidence', 'PMID', 'Source', 'taxonomy'
        regulation = self._build_regulations()
        regulation = self.select_columns(regulation,
                                         'TF_locusTag', 'TG_locusTag', 'Role', 'Evidence', 'PMID', 'Source', 'taxonomy')

        # 'regulator', 'TF_locusTag'
        regulator = self._transform_regulator()

        # 'TF_locusTag', 'TG_locusTag', 'Role', 'Evidence', 'PMID', 'Source', 'regulator', 'taxonomy'
        regulation_regulator = pd.merge(regulation, regulator, on='TF_locusTag')

        # 'operon', 'genes', 'tfbss', 'Operon', 'Orientation', 'Genes'
        operon = self._transform_operon()

        # 'operon', 'Genes'
        operon_by_gene = operon.explode(column='Genes')
        operon_by_gene = self.select_columns(operon_by_gene, 'operon', 'Genes')

        # 'operon', 'Genes'
        # 'TF_locusTag', 'TG_locusTag', 'Role', 'Evidence', 'PMID', 'Source', 'regulator', 'taxonomy'
        regulation_regulator_operon = pd.merge(regulation_regulator, operon_by_gene,
                                               left_on='TG_locusTag', right_on='Genes')
        regulation_regulator_operon = regulation_regulator_operon.drop(columns=['Genes'])

        # 'operon', 'genes', 'tfbss', 'Operon', 'Orientation', 'Genes', 'TF_locusTag', 'TG_locusTag',
        # 'Role', 'Evidence', 'PMID', 'Source', 'taxonomy', 'regulator', 'regulatory_effect'
        regulatory_interaction = pd.merge(regulation_regulator_operon, operon, on='operon')
        regulatory_interaction['regulatory_effect'] = regulatory_interaction['Role']

        regulatory_interaction = apply_processors(regulatory_interaction,
                                                  regulatory_effect=regulatory_effect_coryneregnet)

        regulatory_interaction['regulator_effector'] = None

        regulatory_interaction = self.regulatory_interaction_hash(regulatory_interaction)

        self._stack_transformed_nodes(regulatory_interaction)

        return regulatory_interaction


class RegulatoryInteractionToRegulatorConnector(CoryneRegNetConnector,
                                                source='coryneregnet',
                                                version='0.0.0',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):
    default_from_node = RegulatoryInteraction
    default_to_node = Regulator
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToOperonConnector(CoryneRegNetConnector,
                                             source='coryneregnet',
                                             version='0.0.0',
                                             from_node=RegulatoryInteraction,
                                             to_node=Operon,
                                             register=True):
    default_from_node = RegulatoryInteraction
    default_to_node = Operon
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToGeneConnector(CoryneRegNetConnector,
                                           source='coryneregnet',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):
    default_from_node = RegulatoryInteraction
    default_to_node = Gene
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


class RegulatoryInteractionToTFBSConnector(CoryneRegNetConnector,
                                           source='coryneregnet',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=TFBS,
                                           register=True):
    default_from_node = RegulatoryInteraction
    default_to_node = TFBS
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


class RegulatorToOperonConnector(CoryneRegNetConnector,
                                 source='coryneregnet',
                                 version='0.0.0',
                                 from_node=Regulator,
                                 to_node=Operon,
                                 register=True):
    default_from_node = Regulator
    default_to_node = Operon
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['regulator'].tolist()
        to_identifiers = rin['operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatorToGeneConnector(CoryneRegNetConnector,
                               source='coryneregnet',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):
    default_from_node = Regulator
    default_to_node = Gene
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


class RegulatorToTFBSConnector(CoryneRegNetConnector,
                               source='coryneregnet',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=TFBS,
                               register=True):
    default_from_node = Regulator
    default_to_node = TFBS
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
