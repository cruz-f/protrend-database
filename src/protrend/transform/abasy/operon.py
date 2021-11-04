import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Operon, Gene
from protrend.transform.abasy.base import AbasyTransformer, AbasyConnector
from protrend.transform.abasy.gene import GeneTransformer
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_list, operon_hash
from protrend.utils import SetList


class OperonTransformer(AbasyTransformer,
                        source='abasy',
                        version='0.0.0',
                        node=Operon,
                        order=90,
                        register=True):

    default_transform_stack = {'gene': 'integrated_gene.json'}
    columns = SetList(['protrend_id', 'name', 'genes', 'strand', 'start', 'stop', 'operon_hash',
                       'id', 'source', 'target', 'Effect', 'Evidence', 'taxonomy', 'target_taxonomy'])

    def _transform_networks(self, networks: pd.DataFrame) -> pd.DataFrame:
        networks = networks.dropna(subset=['target', 'taxonomy'])
        networks['target_taxonomy'] = networks['target'] + networks['taxonomy']

        networks = self.drop_duplicates(df=networks, subset=['target_taxonomy'], perfect_match=True)
        networks = networks.dropna(subset=['target_taxonomy'])

        networks = apply_processors(networks,
                                    target_taxonomy=[rstrip, lstrip],
                                    target=[rstrip, lstrip])

        return networks

    def transform(self):
        networks = self._build_networks()
        genes = self._transform_networks(networks)

        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns,
                               reader=read_json_frame)
        gene = self.select_columns(gene, 'protrend_id', 'name', 'strand', 'start', 'stop', 'gene_name_taxonomy')
        gene = gene.rename(columns={'protrend_id': 'genes'})
        gene = apply_processors(gene, genes=to_list)

        df = pd.merge(genes, gene, left_on='target_taxonomy', right_on='gene_name_taxonomy')

        df['tfbss'] = None
        df['operon_hash'] = df['genes']
        df = apply_processors(df, operon_hash=[to_list, operon_hash])
        df = self.drop_duplicates(df=df, subset=['operon_hash'], perfect_match=True)
        df = df.dropna(subset=['operon_hash'])

        self._stack_transformed_nodes(df)

        return df


class OperonToGeneConnector(AbasyConnector,
                            source='abasy',
                            version='0.0.0',
                            from_node=Operon,
                            to_node=Gene,
                            register=False):

    default_connect_stack = {'operon': 'integrated_operon.json'}

    def connect(self):
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)

        operon = apply_processors(operon, genes=to_list)
        operon = operon.explode('genes')

        operon = operon.dropna(subset=['protrend_id', 'genes'])
        operon = operon.drop_duplicates(subset=['protrend_id', 'genes'])

        from_identifiers = operon['protrend_id'].tolist()
        to_identifiers = operon['genes'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
