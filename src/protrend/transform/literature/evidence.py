import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Evidence, Regulator, Operon, Gene
from protrend.transform.literature.base import LiteratureTransformer, LiteratureConnector
from protrend.transform.literature.operon import OperonTransformer
from protrend.transform.literature.regulator import RegulatorTransformer
from protrend.transform.literature.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.utils.processors import apply_processors, to_set_list
from protrend.utils import SetList


class EvidenceTransformer(LiteratureTransformer):
    default_node = Evidence
    default_order = 100
    columns = SetList(['protrend_id',
                       'name', 'description',
                       'regulator_locus_tag', 'operon', 'genes_locus_tag',
                       'regulatory_effect', 'evidence', 'effector', 'mechanism',
                       'publication', 'taxonomy', 'source', 'network_id'])

    def _transform_evidence(self, network: pd.DataFrame) -> pd.DataFrame:
        network = apply_processors(network, evidence=to_set_list)
        network = network.explode(column='evidence')

        network = self.drop_duplicates(df=network, subset=['evidence'], perfect_match=True, preserve_nan=True)
        network = network.dropna(subset=['evidence'])

        def split_evidence(item: str) -> SetList:
            res = SetList()
            comma_split = item.split(',')

            for element in comma_split:
                and_split = element.split(' and ')

                for sub_element in and_split:
                    sub_element = sub_element.rstrip().lstrip()
                    res.append(sub_element)

            return res

        network = apply_processors(network, evidence=split_evidence)
        network = network.explode(column='evidence')

        network = self.drop_duplicates(df=network, subset=['evidence'], perfect_match=True, preserve_nan=True)
        network = network.dropna(subset=['evidence'])

        network['name'] = network['evidence']
        network['description'] = None
        return network

    def transform(self):
        network = self._build_network()
        evidence = self._transform_evidence(network)

        self._stack_transformed_nodes(evidence)
        return evidence


class EvidenceToRegulatorConnector(LiteratureConnector):
    default_from_node = Evidence
    default_to_node = Regulator
    default_connect_stack = {'evidence': 'integrated_evidence.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)

        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)

        df = pd.merge(regulator, evidence, on='network_id', suffixes=('_regulator', '_evidence'))

        from_identifiers = df['protrend_id_evidence'].tolist()
        to_identifiers = df['protrend_id_regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EvidenceToOperonConnector(LiteratureConnector):
    default_from_node = Evidence
    default_to_node = Operon
    default_connect_stack = {'evidence': 'integrated_evidence.json', 'operon': 'integrated_operon.json'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)

        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)

        df = pd.merge(operon, evidence, on='network_id', suffixes=('_operon', '_evidence'))

        from_identifiers = df['protrend_id_evidence'].tolist()
        to_identifiers = df['protrend_id_operon'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EvidenceToGeneConnector(LiteratureConnector):
    default_from_node = Evidence
    default_to_node = Gene
    default_connect_stack = {'evidence': 'integrated_evidence.json', 'operon': 'integrated_operon.json'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)

        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)

        df = pd.merge(operon, evidence, on='network_id', suffixes=('_operon', '_evidence'))
        df = df.explode(column='genes')

        from_identifiers = df['protrend_id_evidence'].tolist()
        to_identifiers = df['genes'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EvidenceToRegulatoryInteractionConnector(LiteratureConnector):
    default_from_node = Evidence
    default_to_node = Gene
    default_connect_stack = {'evidence': 'integrated_evidence.json', 'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        evidence = read_from_stack(stack=self._connect_stack, file='evidence',
                                   default_columns=EvidenceTransformer.columns, reader=read_json_frame)

        rin = read_from_stack(stack=self._connect_stack, file='rin',
                                 default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)

        df = pd.merge(rin, evidence, on='network_id', suffixes=('_rin', '_evidence'))

        from_identifiers = df['protrend_id_evidence'].tolist()
        to_identifiers = df['protrend_id_rin'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
