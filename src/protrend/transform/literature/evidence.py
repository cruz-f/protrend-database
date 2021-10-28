import pandas as pd

from protrend.model.model import Evidence
from protrend.transform.literature.base import LiteratureTransformer
from protrend.transform.processors import apply_processors, to_set_list
from protrend.utils import SetList


class EvidenceTransformer(LiteratureTransformer):
    default_node = Evidence
    default_order = 100
    columns = SetList(['protrend_id',
                       'name', 'description',
                       'regulator_locus_tag', 'regulator_name', 'operon', 'genes_locus_tag',
                       'genes_name', 'regulatory_effect', 'evidence', 'effector', 'mechanism',
                       'publication', 'taxonomy', 'source'])

    def _transform_evidence(self, network: pd.DataFrame) -> pd.DataFrame:
        network = apply_processors(network, evidence=to_set_list)
        network = network.explode(column='evidence')

        network = self.drop_duplicates(df=network, subset=['evidence'], perfect_match=True, preserve_nan=True)
        network = network.dropna(subset=['evidence'])

        def split_evidence(item: str) -> list:
            res = SetList()
            comma_split = item.split(',')

            for element in comma_split:
                and_split = item.split(' and ')

                for sub_element in and_split:
                    sub_element = sub_element.rstrip().lstrip()
                    res.append(sub_element)

            return res

        network = apply_processors(network, evidence=split_evidence)
        network = network.explode(column='evidence')

        network = self.drop_duplicates(df=network, subset=['evidence'], perfect_match=True, preserve_nan=True)
        network = network.dropna(subset=['evidence'])

        network['name'] = network['evidence']
        network['description'] = None
        return network

    def transform(self):
        network = self._build_network()
        evidence = self._transform_evidence(network)

        self._stack_transformed_nodes(evidence)
        return evidence
