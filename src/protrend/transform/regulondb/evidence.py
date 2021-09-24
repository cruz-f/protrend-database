from typing import Union

import pandas as pd

from protrend.io import read_from_stack, read_txt
from protrend.model.model import Evidence
from protrend.transform.processors import apply_processors, rstrip, lstrip
from protrend.transform.regulondb.base import RegulondbTransformer
from protrend.utils import SetList


class EvidenceTransformer(RegulondbTransformer):
    default_node = Evidence
    default_transform_stack = {'evidence': 'evidence.txt'}
    default_order = 100
    columns = SetList(['protrend_id',
                       'name', 'description',
                       'evidence_id', 'evidence_name', 'type_object', 'evidence_code', 'evidence_note',
                       'evidence_internal_comment', 'key_id_org', 'evidence_type', 'evidence_category', 'head',
                       'example'])
    read_columns = SetList(['evidence_id', 'evidence_name', 'type_object', 'evidence_code', 'evidence_note',
                            'evidence_internal_comment', 'key_id_org', 'evidence_type', 'evidence_category', 'head',
                            'example'])

    def _transform_evidence(self, evidence: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=evidence, subset=['evidence_id', 'evidence_name'],
                                  perfect_match=False, preserve_nan=True)

        def remove_evidence_note(item: str) -> Union[None, str]:
            item = item.replace('EVIDENCE_NOTE', '')
            if not item:
                return
            return item

        df = apply_processors(df,
                              evidence_id=[rstrip, lstrip], evidence_name=[rstrip, lstrip],
                              evidence_note=[rstrip, lstrip, remove_evidence_note])
        df = df.dropna(subset=['evidence_id', 'evidence_name'])
        df['name'] = df['evidence_name']
        df['description'] = df['evidence_note']

        return df

    def transform(self):
        evidence = read_from_stack(stack=self.transform_stack, file='evidence', default_columns=self.read_columns,
                                   reader=read_txt, skiprows=38, names=self.read_columns)
        evidence = self._transform_evidence(evidence)

        self._stack_transformed_nodes(evidence)
        return evidence
