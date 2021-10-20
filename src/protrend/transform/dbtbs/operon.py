from statistics import mode, StatisticsError

import numpy as np
import pandas as pd

from protrend.io import read_from_stack, read_json_lines, read_json_frame
from protrend.model.model import Operon, Gene, TFBS
from protrend.transform.dbtbs.base import DBTBSTransformer, DBTBSConnector
from protrend.transform.dbtbs.gene import GeneTransformer
from protrend.transform.dbtbs.tfbs import TFBSTransformer
from protrend.transform.processors import (apply_processors, flatten_set_list, to_list,
                                           to_set_list, operon_hash, to_list_nan)
from protrend.utils import SetList
from protrend.utils.miscellaneous import is_null


class OperonTransformer(DBTBSTransformer):
    default_node = Operon
    default_transform_stack = {'operon': 'Operon.json', 'gene': 'integrated_gene.json', 'tfbs': 'integrated_tfbs.json'}
    default_order = 80
    columns = SetList(['name', 'promoters', 'genes', 'tfbss', 'strand', 'start', 'stop', 'operon_hash', 'protrend_id',
                       'tf', 'url', 'evidence', 'pubmed', 'comment', 'gene', 'tfbs'])
    read_columns = SetList(['name', 'tf', 'url', 'evidence', 'pubmed', 'comment', 'gene', 'tfbs'])

    def _transform_operon(self, operon: pd.DataFrame) -> pd.DataFrame:
        operon = operon.explode(column='name')
        operon = self.drop_duplicates(df=operon, subset=['name'], perfect_match=True, preserve_nan=True)
        operon = operon.dropna(subset=['name'])
        operon = apply_processors(operon, tf=to_list_nan, url=to_list_nan, evidence=to_list_nan, pubmed=to_list_nan,
                                  comment=to_list_nan, gene=to_list_nan, tfbs=to_list_nan)
        return operon

    def _transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = self.select_columns(gene,
                                   'protrend_id', 'locus_tag', 'name', 'strand', 'start', 'stop', 'operon')

        gene = gene.dropna(subset=['protrend_id', 'operon'])
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id',
                                    'locus_tag': 'gene_locus_tag',
                                    'name': 'gene_name',
                                    'operon': 'gene_operon',
                                    'strand': 'gene_strand',
                                    'start': 'gene_start',
                                    'stop': 'gene_stop'})

        return gene

    def _transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = self.select_columns(tfbs, 'protrend_id', 'operon')
        tfbs = tfbs.dropna(subset=['protrend_id', 'operon'])
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id', 'operon': 'tfbs_operon'})
        return tfbs

    def _transform_operon_by_gene(self, operon: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        df = pd.merge(operon, gene, left_on='name', right_on='gene_operon')

        # group by the operon_id
        aggregation = {'gene_protrend_id': to_set_list,
                       'gene_locus_tag': to_set_list,
                       'gene_name': to_set_list,
                       'gene_operon': to_set_list,
                       'gene_strand': to_set_list,
                       'gene_start': to_set_list,
                       'gene_stop': to_set_list}
        df = self.group_by(df=df, column='name', aggregation=aggregation, default=flatten_set_list)

        df = df.rename(columns={'gene_protrend_id': 'genes'})

        df['operon_hash'] = df['genes']

        df = apply_processors(df, operon_hash=[to_list, operon_hash])
        df = df.dropna(subset=['operon_hash'])
        df = self.drop_duplicates(df=df, subset=['operon_hash'], perfect_match=True, preserve_nan=True)

        return df

    def _transform_operon_by_tfbs(self, operon: pd.DataFrame, tfbs: pd.DataFrame) -> pd.DataFrame:

        tfbs_by_operon = tfbs.explode('tfbs_operon')
        df = pd.merge(operon, tfbs_by_operon, how='left', left_on='name', right_on='tfbs_operon')

        aggregation = {'tfbs_protrend_id': to_set_list, 'tfbs_operon': to_set_list}
        df = self.group_by(df=df, column='name', aggregation=aggregation, default=flatten_set_list)

        df = df.rename(columns={'tfbs_protrend_id': 'tfbss'})
        df = df.drop(columns=['tfbs_operon'])

        return df

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
        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=self.read_columns, reader=read_json_lines)
        operon = self._transform_operon(operon)

        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = self._transform_gene(gene)

        tfbs = read_from_stack(stack=self.transform_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = self._transform_tfbs(tfbs)

        # genes
        # tfbss
        # strand
        # start
        # stop
        df = self._transform_operon_by_tfbs(operon=operon, tfbs=tfbs)
        df = self._transform_operon_by_gene(operon=df, gene=gene)
        df = self._operon_coordinates(operon=df)

        self._stack_transformed_nodes(df)

        return df


class OperonToGeneConnector(DBTBSConnector):
    default_from_node = Operon
    default_to_node = Gene
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


class OperonToTFBSConnector(DBTBSConnector):
    default_from_node = Operon
    default_to_node = TFBS
    default_connect_stack = {'operon': 'integrated_operon.json'}

    def connect(self):
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)

        operon = apply_processors(operon, tfbss=to_list)
        operon = operon.explode('tfbss')

        operon = operon.dropna(subset=['protrend_id', 'tfbss'])
        operon = operon.drop_duplicates(subset=['protrend_id', 'tfbss'])

        from_identifiers = operon['protrend_id'].tolist()
        to_identifiers = operon['tfbss'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class GeneToTFBSConnector(DBTBSConnector):
    default_from_node = Gene
    default_to_node = TFBS
    default_connect_stack = {'operon': 'integrated_operon.json'}

    def connect(self):
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)

        operon = apply_processors(operon, genes=to_list, tfbss=to_list)
        operon = operon.explode('tfbss')
        operon = operon.explode('genes')

        operon = operon.dropna(subset=['genes', 'tfbss'])
        operon = operon.drop_duplicates(subset=['genes', 'tfbss'])

        from_identifiers = operon['genes'].tolist()
        to_identifiers = operon['tfbss'].tolist()
        kwargs = dict(operon=operon['protrend_id'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)