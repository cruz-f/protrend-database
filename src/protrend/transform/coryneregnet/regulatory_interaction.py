import pandas as pd

from protrend.io.json import read_json_frame
from protrend.io.utils import read_from_stack
from protrend.model.model import RegulatoryInteraction
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer
from protrend.transform.coryneregnet.operon import OperonTransformer
from protrend.transform.coryneregnet.regulator import RegulatorTransformer
from protrend.transform.processors import (apply_processors, to_list, regulatory_effect_coryneregnet)
from protrend.utils import SetList


class RegulatoryInteractionTransformer(CoryneRegNetTransformer):
    default_node = RegulatoryInteraction
    default_transform_stack = {'regulator': 'integrated_regulator.json',
                               'operon': 'integrated_operon.json'}
    default_order = 70
    columns = SetList(['regulator_effector', 'regulator', 'operon', 'genes', 'tfbss', 'regulatory_effect',
                       'regulatory_interaction_hash', 'protrend_id',
                       'Operon', 'Orientation', 'Genes', 'TF_locusTag', 'TG_locusTag',
                       'Role', 'Evidence', 'PMID', 'Source'])

    def _transform_regulator(self) -> pd.DataFrame:
        regulator = read_from_stack(stack=self.transform_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = self.select_columns(regulator, 'protrend_id', 'TF_locusTag')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def _transform_operon(self) -> pd.DataFrame:
        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = self.select_columns(operon, 'protrend_id', 'genes', 'tfbss', 'Operon')
        operon = operon.rename(columns={'protrend_id': 'operon'})
        operon = apply_processors(operon, genes=to_list, tfbss=to_list)
        return operon

    def transform(self) -> pd.DataFrame:
        # 'TF_locusTag', 'TG_locusTag', 'Role', 'Evidence', 'PMID', 'Source'
        regulation = self._build_regulations()
        regulation = self.select_columns(regulation, 'TF_locusTag', 'TG_locusTag', 'Role', 'Evidence', 'PMID', 'Source')

        # 'regulator', 'TF_locusTag'
        regulator = self._transform_regulator()

        # 'TF_locusTag', 'TG_locusTag', 'Role', 'Evidence', 'PMID', 'Source', 'regulator'
        regulation_regulator = pd.merge(regulation, regulator, on='TF_locusTag')

        # 'Operon', 'Orientation', 'Genes'
        operons = self._build_operons()
        operons = operons.explode(column='Genes')

        # 'Operon', 'Orientation', 'Genes', 'TF_locusTag', 'TG_locusTag', 'Role', 'Evidence', 'PMID', 'Source',
        # 'regulator'
        regulation_regulator_operons = pd.merge(regulation_regulator, operons, left_on='TG_locusTag', right_on='Genes')

        # 'operon', 'genes', 'tfbss', 'Operon'
        operon = self._transform_operon()

        # 'operon', 'genes', 'tfbss', 'Operon', 'Orientation', 'Genes', 'TF_locusTag', 'TG_locusTag',
        # 'Role', 'Evidence', 'PMID', 'Source', 'regulator', 'regulatory_effect'
        regulatory_interaction = pd.merge(regulation_regulator_operons, operon, on='Operon')
        regulatory_interaction['regulatory_effect'] = regulatory_interaction['Role']

        regulatory_interaction = apply_processors(regulatory_interaction,
                                                  regulatory_effect=regulatory_effect_coryneregnet)

        regulatory_interaction = self.regulatory_interaction_hash(regulatory_interaction)

        self._stack_transformed_nodes(regulatory_interaction)

        return regulatory_interaction
