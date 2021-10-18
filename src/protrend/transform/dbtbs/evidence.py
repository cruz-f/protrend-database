import pandas as pd

from protrend.io import read_from_stack, read_json_lines
from protrend.model.model import Evidence
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.transform.processors import apply_processors
from protrend.utils import SetList


class EvidenceTransformer(DBTBSTransformer):
    default_node = Evidence
    default_transform_stack = {'operon': 'Operon.json'}
    default_order = 100
    columns = SetList(['protrend_id', 'name', 'description',
                       'operon', 'tf', 'url', 'evidence', 'pubmed', 'comment', 'gene', 'tfbs'])
    read_columns = SetList(['name', 'tf', 'url', 'evidence', 'pubmed', 'comment', 'gene', 'tfbs'])

    def _transform_evidence(self, evidence: pd.DataFrame) -> pd.DataFrame:
        df = evidence.dropna(subset=['evidence'])

        def split_dbtbs_evidence(item: list) -> list:
            item = item[0]
            return item.split('; ')

        df = apply_processors(df, evidence=split_dbtbs_evidence)
        df = df.explode(column='evidence')

        df = self.drop_duplicates(df=evidence, subset=['evidence'],
                                  perfect_match=True, preserve_nan=True)

        df = df.rename(columns={'name': 'operon'})

        df['name'] = df['evidence']
        df['description'] = None

        return df

    def transform(self):
        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=self.read_columns, reader=read_json_lines)
        evidence = self._transform_evidence(operon)

        self._stack_transformed_nodes(evidence)
        return evidence
