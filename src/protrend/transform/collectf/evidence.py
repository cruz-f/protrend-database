import pandas as pd

from protrend.io import read_from_stack, read_json_lines, read_json_frame
from protrend.model.model import Evidence, Regulator, Operon, TFBS, Gene, RegulatoryInteraction
from protrend.transform.collectf.base import CollectfTransformer, CollectfConnector
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.transform.collectf.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.transform.collectf.tfbs import TFBSTransformer
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_list
from protrend.utils import SetList


class EvidenceTransformer(CollectfTransformer):
    default_node = Evidence
    default_transform_stack = {'evidence': 'ExperimentalEvidence.json'}
    default_order = 100
    columns = SetList(['exp_id', 'regulon', 'tfbs', 'name', 'description', 'protrend_id'])
    read_columns = SetList(['exp_id', 'regulon', 'tfbs'])

    def _transform_evidence(self, evidence: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=evidence, subset=['exp_id'], perfect_match=True, preserve_nan=True)
        df = apply_processors(df, exp_id=[rstrip, lstrip])
        df = df.dropna(subset=['exp_id'])
        df['name'] = df['exp_id']
        df['description'] = None

        return df

    def transform(self):
        evidence = read_from_stack(stack=self.transform_stack, file='evidence',
                                   default_columns=self.read_columns, reader=read_json_lines)
        evidence = self._transform_evidence(evidence)

        self._stack_transformed_nodes(evidence)
        return evidence


class EvidenceToRegulatorConnector(CollectfConnector):
    default_from_node = Evidence
    default_to_node = Regulator
    default_connect_stack = {'evidence': 'integrated_evidence.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)
        evidence = evidence.rename(columns={'protrend_id': 'evidence_protrend_id'})
        evidence = apply_processors(evidence, regulon=to_list)
        evidence = evidence.explode(column='regulon')

        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        df = pd.merge(evidence, regulator, left_on='regulon', right_on='uniprot_accession')
        df = df.drop_duplicates(subset=['evidence_protrend_id', 'regulator_protrend_id'])

        from_identifiers = df['evidence_protrend_id'].tolist()
        to_identifiers = df['regulator_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EvidenceToTFBSConnector(CollectfConnector):
    default_from_node = Evidence
    default_to_node = TFBS
    default_connect_stack = {'evidence': 'integrated_evidence.json',
                             'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)
        evidence = evidence.rename(columns={'protrend_id': 'evidence_protrend_id'})

        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id'})
        tfbs = apply_processors(tfbs, experimental_evidence=to_list)
        tfbs = tfbs.explode(column='experimental_evidence')
        tfbs = apply_processors(tfbs, experimental_evidence=[rstrip, lstrip])

        df = pd.merge(evidence, tfbs, left_on='name', right_on='experimental_evidence')
        df = df.drop_duplicates(subset=['evidence_protrend_id', 'tfbs_protrend_id'])

        from_identifiers = df['evidence_protrend_id'].tolist()
        to_identifiers = df['tfbs_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EvidenceToOperonConnector(CollectfConnector):
    default_from_node = Evidence
    default_to_node = Operon
    default_connect_stack = {'evidence': 'integrated_evidence.json',
                             'rin': 'integrated_regulatoryinteraction.json',
                             'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)
        evidence = evidence.rename(columns={'protrend_id': 'evidence_protrend_id'})

        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        rin = rin.rename(columns={'protrend_id': 'regulatory_interaction_protrend_id',
                                  'operon': 'operon_protrend_id'})
        rin = rin.explode(column='tfbss')

        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id'})
        tfbs = apply_processors(tfbs, experimental_evidence=to_list)
        tfbs = tfbs.explode(column='experimental_evidence')
        tfbs = apply_processors(tfbs, experimental_evidence=[rstrip, lstrip])

        df = pd.merge(evidence, tfbs, left_on='name', right_on='experimental_evidence')
        df = df.drop_duplicates(subset=['evidence_protrend_id', 'tfbs_protrend_id'])

        df = pd.merge(df, rin, left_on='tfbs_protrend_id', right_on='tfbss')
        df = df.drop_duplicates(subset=['evidence_protrend_id', 'operon_protrend_id'])

        from_identifiers = df['evidence_protrend_id'].tolist()
        to_identifiers = df['operon_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EvidenceToGeneConnector(CollectfConnector):
    default_from_node = Evidence
    default_to_node = Gene
    default_connect_stack = {'evidence': 'integrated_evidence.json',
                             'rin': 'integrated_regulatoryinteraction.json',
                             'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)
        evidence = evidence.rename(columns={'protrend_id': 'evidence_protrend_id'})

        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        rin = rin.rename(columns={'protrend_id': 'regulatory_interaction_protrend_id',
                                  'operon': 'operon_protrend_id',
                                  'genes': 'genes_protrend_id'})
        rin = rin.explode(column='tfbss')

        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id'})
        tfbs = apply_processors(tfbs, experimental_evidence=to_list)
        tfbs = tfbs.explode(column='experimental_evidence')
        tfbs = apply_processors(tfbs, experimental_evidence=[rstrip, lstrip])

        df = pd.merge(evidence, tfbs, left_on='name', right_on='experimental_evidence')
        df = df.drop_duplicates(subset=['evidence_protrend_id', 'tfbs_protrend_id'])

        df = pd.merge(df, rin, left_on='tfbs_protrend_id', right_on='tfbss')
        df = df.drop_duplicates(subset=['evidence_protrend_id', 'operon_protrend_id'])
        df = df.explode(column='genes_protrend_id')
        df = df.drop_duplicates(subset=['evidence_protrend_id', 'genes_protrend_id'])

        from_identifiers = df['evidence_protrend_id'].tolist()
        to_identifiers = df['genes_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EvidenceToRegulatoryInteractionConnector(CollectfConnector):
    default_from_node = Evidence
    default_to_node = RegulatoryInteraction
    default_connect_stack = {'evidence': 'integrated_evidence.json',
                             'rin': 'integrated_regulatoryinteraction.json',
                             'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)
        evidence = evidence.rename(columns={'protrend_id': 'evidence_protrend_id'})

        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        rin = rin.rename(columns={'protrend_id': 'regulatory_interaction_protrend_id',
                                  'operon': 'operon_protrend_id',
                                  'genes': 'genes_protrend_id'})
        rin = rin.explode(column='tfbss')

        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id'})
        tfbs = apply_processors(tfbs, experimental_evidence=to_list)
        tfbs = tfbs.explode(column='experimental_evidence')
        tfbs = apply_processors(tfbs, experimental_evidence=[rstrip, lstrip])

        df = pd.merge(evidence, tfbs, left_on='name', right_on='experimental_evidence')
        df = df.drop_duplicates(subset=['evidence_protrend_id', 'tfbs_protrend_id'])

        df = pd.merge(df, rin, left_on='tfbs_protrend_id', right_on='tfbss')
        df = df.drop_duplicates(subset=['evidence_protrend_id', 'regulatory_interaction_protrend_id'])

        from_identifiers = df['evidence_protrend_id'].tolist()
        to_identifiers = df['regulatory_interaction_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)