import pandas as pd

from protrend.io.json import read_json_frame
from protrend.io.utils import read_from_stack
from protrend.model.model import RegulatoryInteraction
from protrend.transform.literature.base import LiteratureTransformer
from protrend.transform.literature.effector import EffectorTransformer
from protrend.transform.literature.operon import OperonTransformer
from protrend.transform.literature.regulator import RegulatorTransformer
from protrend.transform.processors import (apply_processors, to_list, regulatory_effect_literature)
from protrend.utils import SetList


class RegulatoryInteractionTransformer(LiteratureTransformer):
    default_node = RegulatoryInteraction
    default_transform_stack = {'effector': 'integrated_effector.json',
                               'regulator': 'integrated_regulator.json',
                               'operon': 'integrated_operon.json'}
    default_order = 80
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
        regulator = self.select_columns(regulator, 'protrend_id', 'network_id')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def _transform_operon(self) -> pd.DataFrame:
        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = self.select_columns(operon, 'protrend_id', 'genes', 'network_id')
        operon = operon.rename(columns={'protrend_id': 'operon'})
        operon = apply_processors(operon, genes=to_list)
        return operon

    def transform(self) -> pd.DataFrame:
        # 'regulator_locus_tag', 'operon', 'genes_locus_tag',
        # 'regulatory_effect', 'evidence', 'effector', 'mechanism',
        # 'publication', 'taxonomy', 'source', 'network_id'
        network = self._build_network()
        network = network.rename(columns={'operon': 'network_operon'})

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
        network_effector_regulator = pd.merge(network_effector, regulator, on='network_id')

        # 'regulator_locus_tag', 'operon', 'genes_locus_tag',
        # 'regulatory_effect', 'evidence', 'effector', 'mechanism',
        # 'publication', 'taxonomy', 'source', 'network_id', 'regulator_effector', 'regulator', 'operon', 'genes'
        regulatory_interaction = pd.merge(network_effector_regulator, operon, on='network_id')

        regulatory_interaction = apply_processors(regulatory_interaction,
                                                  regulatory_effect=regulatory_effect_literature)

        regulatory_interaction = self.regulatory_interaction_hash(regulatory_interaction)

        self._stack_transformed_nodes(regulatory_interaction)

        return regulatory_interaction
