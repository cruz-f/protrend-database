from typing import Union

import pandas as pd

from protrend.io import read_from_stack
from protrend.model import Evidence, RegulatoryInteraction
from protrend.transform.regulondb.base import RegulondbTransformer, RegulondbConnector, regulondb_reader
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip, split_semi_colon, to_list_nan, remove_html_tags


class EvidenceTransformer(RegulondbTransformer,
                          source='regulondb',
                          version='0.0.0',
                          node=Evidence,
                          order=100,
                          register=True):
    default_transform_stack = {'evidence': 'evidence.txt'}
    columns = SetList(['protrend_id', 'name', 'description',
                       'evidence_id', 'evidence_name', 'type_object', 'evidence_code', 'evidence_note',
                       'evidence_internal_comment', 'key_id_org', 'evidence_type', 'evidence_category', 'head',
                       'example'])
    read_columns = SetList(['evidence_id', 'evidence_name', 'type_object', 'evidence_code', 'evidence_note',
                            'evidence_internal_comment', 'key_id_org', 'evidence_type', 'evidence_category', 'head',
                            'example'])

    def transform_evidence(self, evidence: pd.DataFrame) -> pd.DataFrame:
        evidence = evidence.assign(name=evidence['evidence_name'].copy(), description=evidence['evidence_note'].copy())

        def remove_evidence_note(item: str) -> Union[None, str]:
            item = item.replace('EVIDENCE_NOTE', '')
            if not item:
                return
            return item

        evidence = apply_processors(evidence,
                                    name=[rstrip, lstrip],
                                    description=[rstrip, lstrip, remove_evidence_note, remove_html_tags])

        evidence = evidence.dropna(subset=['evidence_id', 'name'])
        evidence = self.drop_empty_string(evidence, 'evidence_id', 'name')
        evidence = self.drop_duplicates(evidence, subset=['evidence_id', 'name'])

        # dropping Author statement evidences
        mask = evidence['name'] != 'Author statement'
        evidence = evidence[mask]

        return evidence

    def transform(self):
        reader = regulondb_reader(skiprows=38, names=self.read_columns)
        evidence = read_from_stack(stack=self.transform_stack, key='evidence', columns=self.read_columns,
                                   reader=reader)
        evidence = self.transform_evidence(evidence)

        self.stack_transformed_nodes(evidence)
        return evidence


class EvidenceToRegulatoryInteractionConnector(RegulondbConnector,
                                               source='regulondb',
                                               version='0.0.0',
                                               from_node=Evidence,
                                               to_node=RegulatoryInteraction,
                                               register=True):
    default_connect_stack = {'evidence': 'integrated_evidence.json',
                             'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        target_processors = {'evidence': [rstrip, lstrip, split_semi_colon, to_list_nan]}
        source_df, target_df = self.transform_stacks(source='evidence',
                                                     target='rin',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors=target_processors)
        target_df = target_df.explode('evidence')

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on='evidence_name', target_on='evidence')

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)
        self.stack_json(df)
