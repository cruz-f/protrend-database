import pandas as pd

from protrend.io import read_from_stack, read_json_lines
from protrend.model import Evidence, TFBS, RegulatoryInteraction
from protrend.transform.collectf.base import CollectfTransformer, CollectfConnector
from protrend.transform.transformations import drop_empty_string, drop_duplicates
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_list_nan


class EvidenceTransformer(CollectfTransformer,
                          source='collectf',
                          version='0.0.1',
                          node=Evidence,
                          order=100,
                          register=True):
    default_transform_stack = {'evidence': 'ExperimentalEvidence.json'}
    columns = SetList(['protrend_id', 'name', 'description', 'exp_id', 'regulon', 'tfbs'])
    read_columns = SetList(['exp_id', 'regulon', 'tfbs'])

    @staticmethod
    def transform_evidence(evidence: pd.DataFrame) -> pd.DataFrame:
        evidence = apply_processors(evidence, exp_id=[rstrip, lstrip])
        evidence = evidence.dropna(subset=['exp_id'])
        evidence = drop_empty_string(evidence, 'exp_id')
        evidence = drop_duplicates(df=evidence, subset=['exp_id'])
        evidence = evidence.assign(name=evidence['exp_id'].copy(), description=None)
        return evidence

    def transform(self):
        evidence = read_from_stack(stack=self.transform_stack, key='evidence',
                                   columns=self.read_columns, reader=read_json_lines)
        evidence = self.transform_evidence(evidence)

        self.stack_transformed_nodes(evidence)
        return evidence


class EvidenceConnector(CollectfConnector, register=False):
    default_connect_stack = {'evidence': 'integrated_evidence.json',
                             'tfbs': 'integrated_tfbs.json',
                             'rin': 'integrated_regulatoryinteraction.json'}

    def _connect(self, source: str, target: str):
        source_df, target_df = self.transform_stacks(source=source,
                                                     target=target,
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors={'experimental_evidence': [to_list_nan]})
        target_df = target_df.explode('experimental_evidence')
        target_df = apply_processors(target_df, experimental_evidence=[rstrip, lstrip])

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='name', target_on='experimental_evidence')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_json(df)


class EvidenceToTFBSConnector(EvidenceConnector,
                              source='collectf',
                              version='0.0.1',
                              from_node=Evidence,
                              to_node=TFBS,
                              register=True):
    def connect(self):
        return self._connect(source='evidence', target='tfbs')


class EvidenceToRegulatoryInteractionConnector(EvidenceConnector,
                                               source='collectf',
                                               version='0.0.1',
                                               from_node=Evidence,
                                               to_node=RegulatoryInteraction,
                                               register=True):

    def connect(self):
        return self._connect(source='evidence', target='rin')
