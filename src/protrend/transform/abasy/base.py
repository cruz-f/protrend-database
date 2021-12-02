from abc import abstractmethod

import pandas as pd

from protrend.io import read_json
from protrend.transform import Transformer, Connector
from protrend.utils import SetList


def read_abasy_network(file_path: str) -> pd.DataFrame:
    network = read_json(file_path)

    network_edges = network['elements']['edges']

    network_data = [edge['data'] for edge in network_edges]

    network_df = pd.DataFrame(network_data)
    return network_df


class AbasyTransformer(Transformer, source='abasy', version='0.0.0', register=False):

    @abstractmethod
    def transform(self):
        pass

    default_network_stack = {'bsub': 'bsub_network.json',
                             'cglu': 'cglu_network.json',
                             'ecol': 'ecol_network.json',
                             'mtub': 'mtub_network.json',
                             'paer': 'paer_network.json',
                             'saur': 'saur_network.json',
                             'scoe': 'scoe_network.json',
                             'spne': 'spne_network.json',
                             'spyo': 'spyo_network.json'}

    default_gene_stack = {'bsub': 'bsub_genes.tsv',
                          'cglu': 'cglu_genes.tsv',
                          'ecol': 'ecol_genes.tsv',
                          'mtub': 'mtub_genes.tsv',
                          'paer': 'paer_genes.tsv',
                          'saur': 'saur_genes.tsv',
                          'scoe': 'scoe_genes.tsv',
                          'spne': 'spne_genes.tsv',
                          'spyo': 'spyo_genes.tsv'}

    default_network_columns = SetList(['id', 'source', 'target', 'Effect', 'Evidence', 'taxonomy'])

    default_gene_columns = SetList(['Gene_name', 'Locus_tag', 'NCBI_gene_ID', 'Uniprot_ID', 'Synonyms',
                                    'Product_function', 'NDA_component', 'taxonomy'])

    taxa_to_organism_code = {'bsub': '224308',
                             'cglu': '196627',
                             'ecol': '511145',
                             'mtub': '83332',
                             'paer': '208964',
                             'saur': '367830',
                             'scoe': '100226',
                             'spne': '406558',
                             'spyo': '186103'}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self._network_stack = self.build_stack(self.default_network_stack)
        self._gene_stack = self.build_stack(self.default_gene_stack)

    @property
    def network_stack(self):
        return self._network_stack

    @property
    def gene_stack(self):
        return self._gene_stack


class AbasyConnector(Connector, source='abasy', version='0.0.0', register=False):

    @abstractmethod
    def connect(self):
        pass
