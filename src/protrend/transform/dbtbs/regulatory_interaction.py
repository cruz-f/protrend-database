import pandas as pd

from protrend.io.json import read_json_frame
from protrend.io.utils import read_from_stack
from protrend.model.model import RegulatoryInteraction
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.transform.dbtbs.operon import OperonTransformer
from protrend.transform.dbtbs.regulator import RegulatorTransformer
from protrend.transform.dbtbs.tfbs import TFBSTransformer
from protrend.transform.processors import (apply_processors, take_first, regulatory_effect_dbtbs)
from protrend.utils import SetList


class RegulatoryInteractionTransformer(DBTBSTransformer):
    default_node = RegulatoryInteraction
    default_transform_stack = {'regulator': 'integrated_regulator.json',
                               'operon': 'integrated_operon.json',
                               'tfbs': 'integrated_tfbs.json'}
    default_order = 80
    columns = SetList([])

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
        operon = self.select_columns(operon, 'protrend_id', 'name')
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

        regulatory_interaction = regulatory_interaction.drop(columns=['regulator_name_dbtbs', 'operon_name',
                                                                      'tf', 'operon'])
        regulatory_interaction = regulatory_interaction.rename(columns={'regulator_protrend_id': 'regulator',
                                                                        'operon_protrend_id': 'operon',
                                                                        'regulation': 'regulatory_effect'})

        regulatory_interaction = apply_processors(regulatory_interaction,
                                                  regulatory_effect=[take_first, regulatory_effect_dbtbs])
        regulatory_interaction['regulator_effector'] = None

        regulatory_interaction = self.regulatory_interaction_hash(regulatory_interaction)

        self._stack_transformed_nodes(regulatory_interaction)

        return regulatory_interaction