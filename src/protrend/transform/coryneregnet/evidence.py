import pandas as pd

from protrend.model import Evidence, RegulatoryInteraction, TFBS
from protrend.transform.coryneregnet.base import (CoryneRegNetTransformer, CoryneRegNetConnector,
                                                  read_coryneregnet_networks)
from protrend.transform.transformations import drop_empty_string, drop_duplicates
from protrend.utils import SetList


class EvidenceTransformer(CoryneRegNetTransformer,
                          source='coryneregnet',
                          version='0.0.0',
                          node=Evidence,
                          order=100,
                          register=True):
    columns = SetList(['protrend_id', 'name', 'description',
                       'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                       'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                       'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence',
                       'PMID', 'Source', 'taxonomy', 'source'])

    @staticmethod
    def transform_evidence(network: pd.DataFrame) -> pd.DataFrame:
        evidence = network.dropna(subset=['Evidence'])
        evidence = drop_empty_string(evidence, 'Evidence')
        evidence = drop_duplicates(df=evidence, subset=['Evidence'])
        evidence = evidence.assign(name=evidence['Evidence'], description=None)
        return evidence

    def transform(self):
        network = read_coryneregnet_networks(self.source, self.version)
        df = self.transform_evidence(network)

        self.stack_transformed_nodes(df)
        return df


class EvidenceToTFBSConnector(CoryneRegNetConnector,
                              source='coryneregnet',
                              version='0.0.0',
                              from_node=Evidence,
                              to_node=TFBS,
                              register=True):
    default_connect_stack = {'evidence': 'integrated_evidence.json', 'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        df = self.create_connection(source='evidence', target='tfbs',
                                    source_on='Evidence', target_on='Evidence')
        self.stack_connections(df)


class EvidenceToRegulatoryInteractionConnector(CoryneRegNetConnector,
                                               source='coryneregnet',
                                               version='0.0.0',
                                               from_node=Evidence,
                                               to_node=RegulatoryInteraction,
                                               register=True):
    default_connect_stack = {'evidence': 'integrated_evidence.json', 'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='evidence', target='rin',
                                    source_on='Evidence', target_on='Evidence')
        self.stack_connections(df)
