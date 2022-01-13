import pandas as pd

from protrend.model import Evidence, RegulatoryInteraction
from protrend.transform.literature.base import LiteratureTransformer, LiteratureConnector, read_literature_networks
from protrend.transform.transformations import drop_empty_string, drop_duplicates
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, to_list_nan, rstrip, lstrip


class EvidenceTransformer(LiteratureTransformer,
                          source='literature',
                          version='0.0.0',
                          node=Evidence,
                          order=100,
                          register=True):
    columns = SetList(['protrend_id', 'name', 'description',
                       'regulator_locus_tag', 'gene_locus_tag',
                       'regulatory_effect', 'evidence', 'effector_name', 'mechanism',
                       'publication', 'taxonomy', 'source'])

    @staticmethod
    def transform_evidence(network: pd.DataFrame) -> pd.DataFrame:
        network = network.assign(name=network['evidence'].copy(), description=None)

        network = apply_processors(network, name=to_list_nan)
        network = network.explode(column='name')

        network = network.dropna(subset=['name'])
        network = drop_empty_string(network, 'name')
        network = drop_duplicates(df=network, subset=['name'])

        def split_evidence(item: str) -> SetList:
            res = SetList()
            comma_split = item.split(',')

            for element in comma_split:
                and_split = element.split(' and ')

                for sub_element in and_split:
                    sub_element = sub_element.rstrip().lstrip()
                    res.append(sub_element)

            return res

        network = apply_processors(network, name=split_evidence)
        network = network.explode(column='name')

        network = network.dropna(subset=['name'])
        network = drop_empty_string(network, 'name')
        network = drop_duplicates(df=network, subset=['name'])
        return network

    def transform(self):
        network = read_literature_networks(source=self.source, version=self.version)
        evidence = self.transform_evidence(network)

        self.stack_transformed_nodes(evidence)
        return evidence


class EvidenceToRegulatoryInteractionConnector(LiteratureConnector,
                                               source='literature',
                                               version='0.0.0',
                                               from_node=Evidence,
                                               to_node=RegulatoryInteraction,
                                               register=True):
    default_connect_stack = {'evidence': 'integrated_evidence.json', 'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        source_df, target_df = self.transform_stacks(source='evidence',
                                                     target='rin',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_on='evidence',
                                                     target_on='evidence',
                                                     source_processors={'evidence': [to_list_nan]},
                                                     target_processors={'evidence': [to_list_nan]})
        source_df = source_df.explode('evidence')
        source_df = apply_processors(source_df, evidence=[rstrip, lstrip])

        target_df = target_df.explode('evidence')
        target_df = apply_processors(target_df, evidence=[rstrip, lstrip])

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='evidence', target_on='evidence')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_connections(df)
