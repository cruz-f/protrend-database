from typing import List

import pandas as pd

from protrend.io.json import read_json_lines, read_json_frame
from protrend.io.utils import read_from_stack
from protrend.model.model import Effector, Source, Organism, Regulator
from protrend.transform.annotation import annotate_effectors
from protrend.transform.dto import EffectorDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, to_int_str, to_list
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.regprecise.source import SourceTransformer
from protrend.utils import SetList


class EffectorTransformer(RegPreciseTransformer):
    default_node = Effector
    default_transform_stack = {'effector': 'Effector.json'}
    default_order = 100
    columns = SetList(['synonyms', 'mechanism', 'kegg_compounds', 'effector_id', 'url',
                       'regulog', 'name', 'protrend_id'])
    read_columns = SetList(['effector_id', 'name', 'url', 'regulog'])

    def _transform_effector(self, effector: pd.DataFrame):
        effector = effector.dropna(subset=['name'])
        effector = self.drop_duplicates(df=effector, subset=['name'], perfect_match=True, preserve_nan=False)

        effector = apply_processors(effector, effector_id=to_int_str, name=[rstrip, lstrip], regulog=to_int_str)

        effector = self.create_input_value(effector, 'name')

        return effector

    @staticmethod
    def _transform_effectors(names: List[str]):
        dtos = [EffectorDTO(input_value=name) for name in names]
        annotate_effectors(dtos=dtos, names=names)

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        effector = read_from_stack(stack=self.transform_stack, file='effector', default_columns=self.read_columns,
                                   reader=read_json_lines)
        effector = self._transform_effector(effector)

        names = effector['input_value'].tolist()
        effectors = self._transform_effectors(names)

        df = pd.merge(effectors, effector, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        df = df.drop(columns=['input_value'])

        self._stack_transformed_nodes(df)

        return df


class EffectorToSourceConnector(RegPreciseConnector):
    default_from_node = Effector
    default_to_node = Source
    default_connect_stack = {'effector': 'integrated_effector.json', 'source': 'integrated_source.json'}

    def connect(self):
        effector = read_from_stack(stack=self._connect_stack, file='effector',
                                   default_columns=EffectorTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source', default_columns=SourceTransformer.columns,
                                 reader=read_json_frame)

        from_identifiers = effector['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=effector['url'].tolist(),
                      external_identifier=effector['effector_id'].tolist(),
                      key=['effector_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)


class EffectorToOrganismConnector(RegPreciseConnector):
    default_from_node = Effector
    default_to_node = Organism
    default_connect_stack = {'effector': 'integrated_effector.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        effector = read_from_stack(stack=self._connect_stack, file='effector',
                                   default_columns=EffectorTransformer.columns, reader=read_json_frame)
        effector = apply_processors(effector, effector_id=to_int_str)
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = apply_processors(regulator, effector=to_list)
        regulator = regulator.explode('effector')
        regulator = apply_processors(regulator, effector=to_int_str)

        merged = pd.merge(effector, regulator, left_on='effector_id', right_on='effector',
                          suffixes=('_effector', '_regulator'))
        merged = merged.dropna(subset=['protrend_id_effector', 'organism_protrend_id'])
        merged = merged.drop_duplicates(subset=['protrend_id_effector', 'organism_protrend_id'])

        from_identifiers = merged['protrend_id_effector'].tolist()
        to_identifiers = merged['organism_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EffectorToRegulatorConnector(RegPreciseConnector):
    default_from_node = Effector
    default_to_node = Regulator
    default_connect_stack = {'effector': 'integrated_effector.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        effector = read_from_stack(stack=self._connect_stack, file='effector',
                                   default_columns=EffectorTransformer.columns, reader=read_json_frame)
        effector = apply_processors(effector, effector_id=to_int_str)
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = apply_processors(regulator, effector=to_list)
        regulator = regulator.explode('effector')
        regulator = apply_processors(regulator, effector=to_int_str)

        merged = pd.merge(effector, regulator, left_on='effector_id', right_on='effector',
                          suffixes=('_effector', '_regulator'))
        merged = merged.dropna(subset=['protrend_id_effector', 'protrend_id_regulator'])
        merged = merged.drop_duplicates(subset=['protrend_id_effector', 'protrend_id_regulator'])

        from_identifiers = merged['protrend_id_effector'].tolist()
        to_identifiers = merged['protrend_id_regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
