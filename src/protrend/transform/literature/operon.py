from statistics import mode, StatisticsError

import numpy as np
import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Operon
from protrend.transform.literature.base import LiteratureTransformer
from protrend.transform.literature.gene import GeneTransformer
from protrend.transform.processors import (apply_processors, operon_name, to_list, operon_hash, to_str)
from protrend.utils import SetList
from protrend.utils.miscellaneous import is_null


class OperonTransformer(LiteratureTransformer):
    default_node = Operon
    default_transform_stack = {'gene': 'integrated_gene.json'}
    default_order = 90
    columns = SetList(['name', 'genes', 'strand', 'start', 'stop', 'operon_hash', 'protrend_id',
                       'regulator_locus_tag', 'regulator_name', 'operon', 'genes_locus_tag',
                       'genes_name', 'regulatory_effect', 'evidence', 'effector', 'mechanism',
                       'publication', 'taxonomy', 'source', 'network_id'])

    def _transform_operon_by_gene(self, network: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        network = apply_processors(network, operon=to_str)
        network = self.drop_duplicates(df=network, subset=['operon', 'taxonomy'], perfect_match=True, preserve_nan=True)
        network = network.dropna(subset=['operon'])

        operon_gene = pd.merge(gene, network, left_on='gene_operon', right_on='operon')

        operon_gene = operon_gene.rename(columns={'gene_protrend_id': 'genes'})

        operon_gene['operon_hash'] = operon_gene['genes']
        operon_gene['name'] = operon_gene['gene_name']

        operon_gene = apply_processors(operon_gene, operon_hash=[to_list, operon_hash], name=[to_list, operon_name])

        operon_gene = self.drop_duplicates(df=operon_gene, subset=['operon_hash'],
                                           perfect_match=True, preserve_nan=True)
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
        gene = self.select_columns(gene, 'protrend_id', 'locus_tag', 'name', 'strand', 'start', 'stop', 'operon')
        gene = apply_processors(gene, operon=to_str)
        gene = gene.dropna(subset=['protrend_id', 'operon'])
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id',
                                    'locus_tag': 'gene_locus_tag',
                                    'name': 'gene_name',
                                    'operon': 'gene_operon',
                                    'strand': 'gene_strand',
                                    'start': 'gene_start',
                                    'stop': 'gene_stop'})
        return gene

    def transform(self):
        network = self._build_network()
        gene = self._transform_gene()

        operon = self._transform_operon_by_gene(network=network, gene=gene)

        operon = self._operon_coordinates(operon=operon)

        self._stack_transformed_nodes(operon)

        return operon