from typing import Union

import pandas as pd

from protrend.model import Effector
from protrend.transform.literature.base import LiteratureTransformer
from protrend.utils import SetList, is_null
from protrend.utils.processors import apply_processors, to_list_nan


class EffectorTransformer(LiteratureTransformer,
                          source='literature',
                          version='0.0.0',
                          node=Effector,
                          order=100,
                          register=True):
    columns = SetList(['protrend_id', 'name', 'kegg_compounds',
                       'regulator_locus_tag', 'gene_locus_tag',
                       'regulatory_effect', 'evidence', 'effector_name', 'mechanism',
                       'publication', 'taxonomy', 'source'])

    def transform_effector(self, network: pd.DataFrame) -> pd.DataFrame:
        network = network.assign(name=network['effector_name'].copy())

        network = apply_processors(network, name=to_list_nan)
        network = network.explode(column='name')

        network = network.dropna(subset=['name'])
        network = self.drop_empty_string(network, 'name')
        network = self.drop_duplicates(df=network, subset=['name'])

        def filter_map_nan(item: str) -> Union[str, None]:

            if item.lower().rstrip().lstrip() == 'nan':
                return

            if item.lower().rstrip().lstrip() == '@':
                return

            if item.lower().rstrip().lstrip() == 'unk':
                return

            return item

        def filter_map_stress(item: str) -> Union[str, None]:

            if is_null(item):
                return

            return item.replace('stress:', '').replace('stress: ', '').replace(' stress:', '').replace(' stress: ', '')

        def split_effectors(item: str) -> SetList:
            items = item.split('|')

            res = SetList()
            for sub_item in items:

                sub_item = sub_item.rstrip().lstrip()

                sub_items = sub_item.split('/')

                for sub_sub_item in sub_items:
                    sub_sub_item = sub_sub_item.rstrip().lstrip()
                    sub_sub_item = filter_map_nan(sub_sub_item)
                    sub_sub_item = filter_map_stress(sub_sub_item)
                    if not is_null(sub_sub_item):
                        sub_sub_item = sub_sub_item.rstrip().lstrip()
                    res.append(sub_sub_item)
            return res

        network = apply_processors(network, name=split_effectors)
        network = network.explode(column='name')

        network = network.dropna(subset=['name'])
        network = self.drop_empty_string(network, 'name')
        network = self.drop_duplicates(df=network, subset=['name'])

        network = self.create_input_value(network, 'name')
        return network

    def transform(self):
        network = self.read_network()

        effectors = self.transform_effector(network)
        annotated_effectors = self.annotate_effectors(effectors)

        df = pd.merge(annotated_effectors, effectors, on='input_value', suffixes=('_annotation', '_literature'))

        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_literature')

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
