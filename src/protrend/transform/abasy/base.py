from abc import abstractmethod

import pandas as pd

from protrend.io import read_json, read_csv
from protrend.transform import MultiStackTransformer, Connector
from protrend.utils import SetList, MultiStack


def read_abasy_network(file_path: str) -> pd.DataFrame:
    network = read_json(file_path)

    network_edges = network['elements']['edges']

    network_data = [edge['data'] for edge in network_edges]

    network_df = pd.DataFrame(network_data)
    network_df = network_df.rename(columns={'source': 'regulator'})
    return network_df


class AbasyTransformer(MultiStackTransformer, source='abasy', version='0.0.0', register=False):
    _taxa = ['224308',
             '196627',
             '511145',
             '83332',
             '208964',
             '367830',
             '100226',
             '406558',
             '186103']
    _source = ['abasy'] * 9
    _net_reader = [read_abasy_network] * 9
    _gene_reader = [read_csv] * 9

    default_transform_stack = {
        'network': MultiStack(
            stack=['bsub_network.json',
                   'cglu_network.json',
                   'ecol_network.json',
                   'mtub_network.json',
                   'paer_network.json',
                   'saur_network.json',
                   'scoe_network.json',
                   'spne_network.json',
                   'spyo_network.json'],
            taxa=_taxa,
            source=_source,
            reader=_net_reader
        ),
        'gene': MultiStack(
            stack=['bsub_genes.tsv',
                   'cglu_genes.tsv',
                   'ecol_genes.tsv',
                   'mtub_genes.tsv',
                   'paer_genes.tsv',
                   'saur_genes.tsv',
                   'scoe_genes.tsv',
                   'spne_genes.tsv',
                   'spyo_genes.tsv'],
            taxa=_taxa,
            source=_source,
            reader=_gene_reader
        )
    }

    default_network_columns = SetList(['id', 'regulator', 'target', 'Effect', 'Evidence', 'source', 'taxonomy'])

    default_gene_columns = SetList(['Gene_name', 'Locus_tag', 'NCBI_gene_ID', 'Uniprot_ID', 'Synonyms',
                                    'Product_function', 'NDA_component', 'taxonomy'])

    @abstractmethod
    def transform(self):
        pass


class AbasyConnector(Connector, source='abasy', version='0.0.0', register=False):

    @abstractmethod
    def connect(self):
        pass
