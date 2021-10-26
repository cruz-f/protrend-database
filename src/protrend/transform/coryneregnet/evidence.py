import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Evidence, RegulatoryInteraction
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer, CoryneRegNetConnector
from protrend.transform.coryneregnet.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.utils import SetList


class EvidenceTransformer(CoryneRegNetTransformer):
    default_node = Evidence
    default_transform_stack = {'bsub': 'bsub_regulation.csv',
                               'cglu': 'cglu_regulation.csv',
                               'cglu_rna': 'cglu_rna.csv',
                               'ecol': 'ecol_regulation.csv',
                               'mtub': 'mtub_regulation.csv'}
    default_order = 100
    columns = SetList(['protrend_id',
                       'name', 'description',
                       'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                       'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                       'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence', 'PMID', 'Source', 'taxonomy'])

    def _transform_evidence(self, regulation: pd.DataFrame) -> pd.DataFrame:
        regulation = self.drop_duplicates(df=regulation, subset=['Evidence'], perfect_match=True, preserve_nan=True)
        regulation = regulation.dropna(subset=['Evidence'])
        regulation['name'] = regulation['Evidence']
        regulation['description'] = None
        return regulation

    def transform(self):
        regulation = self._build_regulations()
        evidence = self._transform_evidence(regulation)

        self._stack_transformed_nodes(evidence)
        return evidence


class EvidenceToRegulatoryInteractionConnector(CoryneRegNetConnector):
    default_from_node = Evidence
    default_to_node = RegulatoryInteraction
    default_connect_stack = {'evidence': 'integrated_evidence.json', 'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)

        df = pd.merge(rin, evidence, on='Evidence', suffixes=('_rin', '_evidence'))

        from_identifiers = df['protrend_id_evidence'].tolist()
        to_identifiers = df['protrend_id_rin'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
