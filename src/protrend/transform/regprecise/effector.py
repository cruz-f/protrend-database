from typing import List

import pandas as pd

from protrend.io.utils import read_from_stack
from protrend.transform.annotation import annotate_effectors
from protrend.transform.connector import DefaultConnector
from protrend.transform.transformer import Transformer
from protrend.transform.dto import EffectorDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors
from protrend.transform.regprecise import RegulatorTransformer, SourceTransformer
from protrend.transform.regprecise.settings import (EffectorSettings, EffectorToSource, EffectorToOrganism,
                                                    EffectorToRegulator)


class EffectorTransformer(Transformer):
    default_settings = EffectorSettings
    columns = {'protrend_id', 'effector_id', 'name', 'url', 'regulog', 'mechanism', 'synonyms', 'kegg_compounds'}
    read_columns = {'effector_id', 'name', 'url', 'regulog'}

    def _transform_effector(self, effector: pd.DataFrame):
        effector = effector.dropna(subset=['name'])
        effector = self.drop_duplicates(df=effector, subset=['name'], perfect_match=True, preserve_nan=False)

        apply_processors(rstrip, lstrip, df=effector, col='name')

        effector = self.create_input_value(effector, 'name')

        return effector

    @staticmethod
    def _transform_effectors(names: List[str]):
        dtos = [EffectorDTO(input_value=name) for name in names]
        annotate_effectors(dtos=dtos, names=names)

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        effector = read_from_stack(tl=self, file='effector', json=True, default_columns=self.read_columns)
        effector = self._transform_effector(effector)

        names = effector['input_value'].tolist()
        effectors = self._transform_effectors(names)

        df = pd.merge(effectors, effector, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        df = df.drop(columns=['input_value'])

        self._stack_transformed_nodes(df)

        return df


class EffectorToSourceConnector(DefaultConnector):
    default_settings = EffectorToSource

    def connect(self):
        effector = read_from_stack(tl=self, file='effector', json=False, default_columns=EffectorTransformer.columns)
        source = read_from_stack(tl=self, file='source', json=False, default_columns=SourceTransformer.columns)

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

        self.stack_csv(df)


class EffectorToOrganismConnector(DefaultConnector):
    default_settings = EffectorToOrganism

    def connect(self):
        effector = read_from_stack(tl=self, file='effector', json=False, default_columns=EffectorTransformer.columns)
        regulator = read_from_stack(tl=self, file='regulator', json=False, default_columns=RegulatorTransformer.columns)
        regulator = regulator.explode('effector')

        merged = pd.merge(effector, regulator, left_on='effector_id', right_on='effector',
                          suffixes=('_effector', '_regulator'))
        merged = merged.dropna(subset=['protrend_id_effector'])
        merged = merged.dropna(subset=['organism_protrend_id'])
        merged = merged.drop_duplicates(subset=['protrend_id_effector', 'organism_protrend_id'])

        from_identifiers = merged['protrend_id_effector'].tolist()
        to_identifiers = merged['organism_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)


class EffectorToRegulatorConnector(DefaultConnector):
    default_settings = EffectorToRegulator

    def connect(self):
        effector = read_from_stack(tl=self, file='effector', json=False, default_columns=EffectorTransformer.columns)
        regulator = read_from_stack(tl=self, file='regulator', json=False, default_columns=RegulatorTransformer.columns)
        regulator = regulator.explode('effector')

        merged = pd.merge(effector, regulator, left_on='effector_id', right_on='effector',
                          suffixes=('_effector', '_regulator'))
        merged = merged.dropna(subset=['protrend_id_effector'])
        merged = merged.dropna(subset=['protrend_id_regulator'])
        merged = merged.drop_duplicates(subset=['protrend_id_effector', 'protrend_id_regulator'])

        from_identifiers = merged['protrend_id_effector'].tolist()
        to_identifiers = merged['protrend_id_regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)
