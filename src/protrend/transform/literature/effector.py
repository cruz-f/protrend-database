from typing import List, Union

import pandas as pd

from protrend.model.model import Effector
from protrend.transform.annotation import annotate_effectors
from protrend.transform.dto import EffectorDTO
from protrend.transform.literature.base import LiteratureTransformer
from protrend.transform.processors import apply_processors, to_set_list
from protrend.utils import SetList
from protrend.utils.miscellaneous import is_null


class EffectorTransformer(LiteratureTransformer):
    default_node = Effector
    default_order = 100
    columns = SetList(['name', 'mechanism', 'kegg_compounds', 'protrend_id'])

    def _transform_effector(self, network: pd.DataFrame) -> pd.DataFrame:
        network = apply_processors(network, effector=to_set_list)
        network = network.explode(column='effector')

        network = self.drop_duplicates(df=network, subset=['effector'], perfect_match=True, preserve_nan=True)
        network = network.dropna(subset=['effector'])

        def _filter_map_nan(item: str) -> Union[str, None]:

            if item.lower().rstrip().lstrip() == 'nan':
                return None

            if item.lower().rstrip().lstrip() == '@':
                return None

            if item.lower().rstrip().lstrip() == 'unk':
                return None

            return item

        def _filter_map_stress(item: str) -> Union[str, None]:

            if is_null(item):
                return None

            return item.replace('stress:', '').replace('stress: ', '').replace(' stress:', '').replace(' stress: ', '')

        def split_effectors(item: str) -> list:
            items = item.split('|')

            res = SetList()
            for item in items:

                sub_items = item.split('/')

                for sub_item in sub_items:

                    sub_item = sub_item.rstrip().lstrip()
                    sub_item = _filter_map_nan(sub_item)
                    sub_item = _filter_map_stress(sub_item)
                    res.append(sub_item)
            return res

        network = apply_processors(network, effector=split_effectors)
        network = network.explode(column='effector')

        network = self.drop_duplicates(df=network, subset=['effector'], perfect_match=True, preserve_nan=True)
        network = network.dropna(subset=['effector'])

        network['name'] = network['effector']
        network['mechanism'] = None

        network = self.create_input_value(network, 'name')

        return network

    @staticmethod
    def _transform_effectors(names: List[str]):
        dtos = [EffectorDTO(input_value=name) for name in names]
        annotate_effectors(dtos=dtos, names=names)

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        network = self._build_network()
        effector = self._transform_effector(network)

        names = effector['input_value'].tolist()
        effectors = self._transform_effectors(names)

        df = pd.merge(effectors, effector, on='input_value', suffixes=('_annotation', '_literature'))

        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_literature')

        df = df.drop(columns=['input_value'])

        self._stack_transformed_nodes(df)

        return df