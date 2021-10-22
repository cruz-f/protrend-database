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
                       ])

    def _transform_operon_by_gene(self, operon: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        # 'Operon', 'Orientation', 'Genes', 'tfbss', 'tfbs_operon'
        operon = operon.explode(column='Genes')

        # 'Operon', 'Orientation', 'Genes', 'tfbss', 'tfbs_operon', 'gene_protrend_id', 'gene_locus_tag', 'gene_name',
        # 'gene_TG_locusTag', 'gene_strand', 'gene_start', 'gene_stop'
        operon_gene = pd.merge(operon, gene, left_on='Genes', right_on='gene_TG_locusTag')

        # group by the operon_id
        aggregation = {'Orientation': take_first, 'tfbss': flatten_set_list}
        operon_gene = self.group_by(df=operon_gene, column='Operon', aggregation=aggregation, default=to_set_list)

        operon_gene = operon_gene.rename(columns={'gene_protrend_id': 'genes'})

        operon_gene['operon_hash'] = operon_gene['genes']
        operon_gene['name'] = operon_gene['gene_name']

        operon_gene = apply_processors(operon_gene, operon_hash=[to_list, operon_hash], name=[to_list, operon_name])

        operon_gene = self.drop_duplicates(df=operon_gene, subset=['operon_hash'],
                                           perfect_match=True, preserve_nan=True)
        operon_gene = operon_gene.dropna(subset=['operon_hash'])

        return operon_gene

    def _transform_operon_by_tfbs(self, operon: pd.DataFrame, tfbs: pd.DataFrame) -> pd.DataFrame:
        # 'Operon', 'Orientation', 'Genes', 'tfbs_protrend_id', 'tfbs_operon'
        operon_tfbs = pd.merge(operon, tfbs, left_on='Operon', right_on='tfbs_operon')

        aggregation = {'Orientation': take_first, 'tfbs_operon': take_first, 'Genes': flatten_set_list}
        operon_tfbs = self.group_by(df=operon_tfbs, column='Operon', aggregation=aggregation, default=to_set_list)

        operon_tfbs = operon_tfbs.rename(columns={'tfbs_protrend_id': 'tfbss'})
        return operon_tfbs

    @staticmethod
    def _operon_coordinates(operon: pd.DataFrame) -> pd.DataFrame:

        def strand(item):
            if is_null(item):
                return None

            if item == '-':
                return 'reverse'

            if item == '+':
                return 'forward'

            return None

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

        operon['strand'] = operon['Orientation'].map(strand, na_action='ignore')
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
        gene = self.select_columns(gene, 'protrend_id', 'locus_tag', 'name', 'strand', 'start', 'stop', 'TG_locusTag')
        gene = gene.dropna(subset=['protrend_id', 'TG_locusTag'])
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id',
                                    'locus_tag': 'gene_locus_tag',
                                    'name': 'gene_name',
                                    'TG_locusTag': 'gene_TG_locusTag',
                                    'strand': 'gene_strand',
                                    'start': 'gene_start',
                                    'stop': 'gene_stop'})
        return gene

    def _transform_tfbs(self) -> pd.DataFrame:
        tfbs = read_from_stack(stack=self.transform_stack, file='tfbs', default_columns=TFBSTransformer.columns,
                               reader=read_json_frame)
        tfbs = self.select_columns(tfbs, 'protrend_id', 'Operon')
        tfbs = tfbs.dropna(subset=['protrend_id', 'Operon'])
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id',
                                    'Operon': 'tfbs_operon'})
        return tfbs

    def transform(self):
        # 'Operon', 'Orientation', 'Genes'
        operon = self._build_operons()

        # 'tfbs_protrend_id', 'tfbs_operon'
        tfbs = self._transform_tfbs()

        # 'gene_protrend_id', 'gene_locus_tag', 'gene_name', 'gene_TG_locusTag', 'gene_strand', 'gene_start',
        # 'gene_stop'
        gene = self._transform_gene()

        # 'Operon', 'Orientation', 'Genes', 'tfbss', 'tfbs_operon'
        operon_tfbs = self._transform_operon_by_tfbs(operon=operon, tfbs=tfbs)
        operon_tfbs_gene = self._transform_operon_by_gene(operon=operon_tfbs, gene=gene)

        operon_tfbs_gene = self._operon_coordinates(operon=operon_tfbs_gene)

        self._stack_transformed_nodes(operon_tfbs_gene)

        return operon_tfbs_gene
