import pandas as pd

from protrend.io import read_from_stack, read_json_lines
from protrend.model.model import Evidence
from protrend.transform.collectf.base import CollectfTransformer
from protrend.transform.processors import apply_processors, remove_white_space, rstrip, lstrip, \
    remove_multiple_white_space, parse_collectf_description


class EvidenceTransformer(CollectfTransformer):
    default_node = Evidence
    default_node_factors = ('name', )
    default_transform_stack = {'evidence': 'ExperimentalEvidence.json'}
    default_order = 100
    columns = {'protrend_id',
               'name', 'description',
               'exp_id', 'regulon', 'tfbs'}
    read_columns = {'exp_id', 'regulon', 'tfbs'}

    def _transform_evidence(self, evidence: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=evidence, subset=['exp_id'], perfect_match=True, preserve_nan=True)

        df = apply_processors(df, exp_id=[rstrip, lstrip])

        df['name'] = df['exp_id']
        df['description'] = None

        return df

    def transform(self):
        evidence = read_from_stack(stack=self.transform_stack, file='evidence',
                                   default_columns=self.read_columns, reader=read_json_lines)
        evidence = self._transform_evidence(evidence)

        self._stack_transformed_nodes(evidence)
        return evidence
