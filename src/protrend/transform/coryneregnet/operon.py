from statistics import mode, StatisticsError

import numpy as np
import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Operon
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer
from protrend.transform.coryneregnet.gene import GeneTransformer
from protrend.transform.coryneregnet.tfbs import TFBSTransformer
from protrend.transform.processors import (apply_processors, operon_name, flatten_set_list, to_list,
                                           to_set_list, operon_hash, take_first)
from protrend.utils import SetList
from protrend.utils.miscellaneous import is_null


class OperonTransformer(CoryneRegNetTransformer):
    default_node = Operon
    default_transform_stack = {'bsub': 'bsub_regulation.csv',
                               'cglu': 'cglu_regulation.csv',
                               'ecol': 'ecol_regulation.csv',
                               'mtub': 'mtub_regulation.csv',
                               'gene': 'integrated_gene.json'}
    default_order = 80
    columns = SetList(['name', 'promoters', 'genes', 'tfbss', 'strand', 'start', 'stop', 'operon_hash', 'protrend_id',
                       'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                       'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                       'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence', 'PMID', 'Source', 'taxonomy',
                       'Orientation', 'Genes'])

    def _transform_operon_by_gene(self, operon: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        operon = operon.explode(column='Genes')

        # 'protrend_id', 'locus_tag', 'name', 'strand', 'start', 'stop',
        # 'TG_locusTag'
        # 'Operon', 'Orientation', 'Genes'
        operon = pd.merge(operon, gene, left_on='Genes', right_on='gene_TG_locusTag')

        # group by the operon_id
        aggregation = {'Orientation': take_first, 'tfbss': flatten_set_list}
        operon = self.group_by(df=operon, column='Operon', aggregation=aggregation, default=to_set_list)

        operon = operon.rename(columns={'gene_protrend_id': 'genes'})

        operon['operon_hash'] = operon['genes']
        operon['name'] = operon['gene_name']

        operon = apply_processors(operon, operon_hash=[to_list, operon_hash], name=[to_list, operon_name])
        operon = operon.dropna(subset=['operon_hash'])
        operon = self.drop_duplicates(df=operon, subset=['operon_hash'], perfect_match=True, preserve_nan=True)

        return operon

    def _transform_operon_by_tfbs(self, operon: pd.DataFrame, tfbs: pd.DataFrame) -> pd.DataFrame:
        operon_tfbs = pd.merge(operon, tfbs, on='Operon')

        aggregation = {'Orientation': take_first, 'Genes': flatten_set_list}
        operon = self.group_by(df=operon, column='Operon', aggregation=aggregation, default=to_set_list)

        operon = operon.rename(columns={'tfbs_protrend_id': 'tfbss'})

        return operon

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

    def transform(self):
        # 'Operon', 'Orientation', 'Genes'
        operon = self._build_operons()

        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = self.select_columns(gene, 'protrend_id', 'locus_tag', 'name', 'strand', 'start', 'stop',
                                   'TG_locusTag')

        gene = gene.dropna(subset=['protrend_id', 'TG_locusTag'])
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id',
                                    'locus_tag': 'gene_locus_tag',
                                    'name': 'gene_name',
                                    'TG_locusTag': 'gene_TG_locusTag',
                                    'strand': 'gene_strand',
                                    'start': 'gene_start',
                                    'stop': 'gene_stop'})

        tfbs = read_from_stack(stack=self.transform_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = self.select_columns(tfbs, 'protrend_id', 'Operon')

        tfbs = tfbs.dropna(subset=['protrend_id'])
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id'})

        # genes
        # promoters
        # tfbss
        # strand
        # start
        # stop
        operon = self._transform_operon_by_tfbs(operon=operon, tfbs=tfbs)
        operon = self._transform_operon_by_gene(operon=operon, gene=gene)
        operon = self._operon_coordinates(operon=operon)

        self._stack_transformed_nodes(operon)

        return operon
