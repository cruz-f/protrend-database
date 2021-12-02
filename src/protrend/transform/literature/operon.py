from statistics import mode, StatisticsError

import numpy as np
import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model import Operon, Gene
from protrend.transform.literature.base import LiteratureTransformer, LiteratureConnector
from protrend.transform.literature.gene import GeneTransformer
from protrend.utils.processors import (apply_processors, operon_name, to_list, operon_hash, to_set_list,
                                       flatten_set_list, take_first)
from protrend.utils import SetList, is_null


class OperonTransformer(LiteratureTransformer,
                        source='literature',
                        version='0.0.0',
                        node=Operon,
                        order=90,
                        register=True):
    default_transform_stack = {'gene': 'integrated_gene.json'}
    columns = SetList(['name', 'genes', 'strand', 'start', 'stop', 'operon_hash', 'protrend_id',
                       'regulator_locus_tag', 'operon', 'genes_locus_tag', 'regulatory_effect', 'evidence',
                       'effector', 'mechanism', 'publication', 'taxonomy', 'source', 'operon_id', 'network_id'])

    def _transform_operon_by_gene(self, network: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        network = self.drop_duplicates(df=network, subset=['operon', 'taxonomy'], perfect_match=True)
        network = network.dropna(subset=['operon'])

        network['operon_id'] = network['operon'] + network['taxonomy']

        operon_gene = pd.merge(gene, network, on='network_id')

        aggregation = {'genes_locus_tag': flatten_set_list, 'operon': take_first, 'taxonomy': take_first,
                       'source': take_first, 'network_id': take_first}
        operon_gene = self.group_by(df=operon_gene, column='operon_id', aggregation=aggregation, default=to_set_list)

        operon_gene = operon_gene.rename(columns={'gene_protrend_id': 'genes'})

        operon_gene['operon_hash'] = operon_gene['genes']
        operon_gene['name'] = operon_gene['gene_name']

        operon_gene = apply_processors(operon_gene, operon_hash=[to_list, operon_hash], name=[to_list, operon_name])

        operon_gene = self.drop_duplicates(df=operon_gene, subset=['operon_hash'], perfect_match=True)
        operon_gene = operon_gene.dropna(subset=['operon_hash'])

        return operon_gene

    @staticmethod
    def _operon_coordinates(operon: pd.DataFrame) -> pd.DataFrame:

        def strand_mode(item):

            if is_null(item):
                return None

            try:
                m = mode(item)

                if is_null(m):
                    return None

                return m

            except StatisticsError:
                for sub_item in item:
                    return sub_item

        def start(item):
            if is_null(item):
                return None

            item = to_list(item)

            x = np.array(item, dtype=np.float64)
            return np.nanmin(x)

        def stop(item):
            if is_null(item):
                return None

            item = to_list(item)

            x = np.array(item, dtype=np.float64)
            return np.nanmax(x)

        operon['strand'] = operon['gene_strand'].map(strand_mode, na_action='ignore')
        forward = operon['strand'] == 'forward'
        reverse = operon['strand'] == 'reverse'

        operon['start'] = None
        operon['stop'] = None

        operon.loc[forward, 'start'] = operon.loc[forward, 'gene_start'].map(start, na_action='ignore')
        operon.loc[forward, 'stop'] = operon.loc[forward, 'gene_stop'].map(stop, na_action='ignore')

        operon.loc[reverse, 'start'] = operon.loc[reverse, 'gene_start'].map(stop, na_action='ignore')
        operon.loc[reverse, 'stop'] = operon.loc[reverse, 'gene_stop'].map(start, na_action='ignore')

        strand_mask = (operon['strand'] != 'reverse') & (operon['strand'] != 'forward')
        operon.loc[strand_mask, 'strand'] = None

        return operon

    def _transform_gene(self) -> pd.DataFrame:
        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = self.select_columns(gene, 'protrend_id', 'locus_tag', 'name', 'strand', 'start', 'stop', 'network_id')
        gene = gene.dropna(subset=['protrend_id', 'network_id'])
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id',
                                    'locus_tag': 'gene_locus_tag',
                                    'name': 'gene_name',
                                    'strand': 'gene_strand',
                                    'start': 'gene_start',
                                    'stop': 'gene_stop'})
        return gene

    def transform(self):
        network = self._build_network()
        network = self.select_columns(network, 'operon', 'genes_locus_tag', 'taxonomy', 'source', 'network_id')
        gene = self._transform_gene()

        operon = self._transform_operon_by_gene(network=network, gene=gene)

        operon = self._operon_coordinates(operon=operon)

        self.stack_transformed_nodes(operon)

        return operon


class OperonToGeneConnector(LiteratureConnector,
                            source='literature',
                            version='0.0.0',
                            from_node=Operon,
                            to_node=Gene,
                            register=True):
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
