import pandas as pd

from protrend.io import read_from_stack, read_json_lines, read_json_frame
from protrend.model.model import Evidence, Operon
from protrend.transform.dbtbs.base import DBTBSTransformer, DBTBSConnector
from protrend.transform.dbtbs.operon import OperonTransformer
from protrend.transform.processors import apply_processors
from protrend.utils import SetList


class EvidenceTransformer(DBTBSTransformer):
    default_node = Evidence
    default_transform_stack = {'operon': 'Operon.json'}
    default_order = 100
    columns = SetList(['protrend_id', 'name', 'description',
                       'operon', 'tf', 'url', 'evidence', 'pubmed', 'comment', 'gene', 'tfbs'])
    read_columns = SetList(['name', 'tf', 'url', 'evidence', 'pubmed', 'comment', 'gene', 'tfbs'])

    def _transform_evidence(self, operon: pd.DataFrame) -> pd.DataFrame:
        df = operon.explode(column='name')

        df = df.dropna(subset=['evidence'])

        def split_dbtbs_evidence(item: list) -> list:
            item = item[0]
            return item.split('; ')

        df = apply_processors(df, evidence=split_dbtbs_evidence)
        df = df.explode(column='evidence')

        df = self.drop_duplicates(df=df, subset=['evidence'], perfect_match=True, preserve_nan=True)

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


class EvidenceToOperonConnector(DBTBSConnector):
    default_from_node = Evidence
    default_to_node = Operon
    default_connect_stack = {'operon': 'integrated_operon.json', 'evidence': 'integrated_evidence.json'}

    def connect(self):
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)

        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)

        df = pd.merge(evidence, operon, left_on='operon', right_on='name', suffixes=('_evidence', '_operon'))

        from_identifiers = df['protrend_id_evidence'].tolist()
        to_identifiers = df['protrend_id_operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
