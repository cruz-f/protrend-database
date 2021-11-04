from statistics import mode, StatisticsError

import numpy as np
import pandas as pd

from protrend.io.json import read_json_lines, read_json_frame
from protrend.io.utils import read_from_stack
from protrend.model.model import Operon, Gene, TFBS
from protrend.transform.collectf.base import CollectfTransformer, CollectfConnector
from protrend.transform.collectf.gene import GeneTransformer
from protrend.transform.collectf.tfbs import TFBSTransformer
from protrend.transform.processors import (apply_processors, operon_name, flatten_set_list, to_list,
                                           to_set_list, operon_hash)
from protrend.utils import build_graph, find_connected_nodes, SetList
from protrend.utils.miscellaneous import is_null


class OperonTransformer(CollectfTransformer):
    default_node = Operon
    default_transform_stack = {'operon': 'Operon.json', 'gene': 'integrated_gene.json', 'tfbs': 'integrated_tfbs.json'}
    default_order = 60
    columns = SetList(['operon_id_new', 'gene', 'operon_id_old', 'regulon', 'tfbs', 'tfbss',
                       'genes', 'gene_locus_tag', 'gene_name', 'gene_old_locus_tag',
                       'gene_strand', 'gene_start', 'gene_stop', 'organism_protrend_id',
                       'operon_id', 'name', 'strand', 'start', 'stop', 'operon_hash', 'protrend_id'])
    read_columns = SetList(['operon_id', 'regulon', 'gene', 'tfbs'])

    def _operon_by_gene(self, operon: pd.DataFrame) -> pd.DataFrame:
        # group duplicates
        operon = operon.explode('gene')

        aggregation = {'operon_id': to_set_list}
        operon = self.group_by(df=operon, column='gene', aggregation=aggregation, default=flatten_set_list)

        return operon

    @staticmethod
    def _normalize_operon(genes) -> pd.DataFrame:

        graph = build_graph(genes)
        operons = find_connected_nodes(graph)
        operon = pd.DataFrame({'gene': operons})

        operon['operon_id'] = operon['gene']
        operon = apply_processors(operon, operon_id=[to_list, operon_hash], gene=to_list)
        operon = operon.explode('gene')
        return operon

    def _transform_operon_by_gene(self, operon: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        operon = apply_processors(operon, gene=to_list)
        genes = operon['gene'].tolist()
        normalized_operon = self._normalize_operon(genes)

        operon_by_gene = self._operon_by_gene(operon)

        operon = pd.merge(normalized_operon, operon_by_gene, on='gene', suffixes=("_new", "_old"))

        # group by the new operons
        aggregation = {'gene': to_set_list}
        operon = self.group_by(df=operon, column='operon_id_new', aggregation=aggregation, default=flatten_set_list)

        operon = apply_processors(operon, gene=to_list)
        operon = operon.explode(column='gene')

        operon = pd.merge(operon, gene, left_on='gene', right_on='gene_old_locus_tag')

        # group by the new genes
        aggregation = {'gene': to_set_list,
                       'gene_protrend_id': to_set_list,
                       'gene_locus_tag': to_set_list,
                       'gene_name': to_set_list,
                       'gene_old_locus_tag': to_set_list,
                       'gene_strand': to_set_list,
                       'gene_start': to_set_list,
                       'gene_stop': to_set_list,
                       'organism_protrend_id': to_set_list}
        operon = self.group_by(df=operon, column='operon_id_new', aggregation=aggregation, default=flatten_set_list)

        mask = operon['organism_protrend_id'].map(len) == 1
        operon = operon[mask]

        operon = operon.rename(columns={'gene_protrend_id': 'genes'})

        operon['operon_hash'] = operon['genes']
        operon['name'] = operon['gene_name']

        operon = apply_processors(operon, operon_hash=[to_list, operon_hash], name=[to_list, operon_name])
        operon = operon.dropna(subset=['operon_hash'])
        operon = self.drop_duplicates(df=operon, subset=['operon_hash'], perfect_match=True, preserve_nan=True)

        return operon

    def _transform_operon_by_tfbs(self, operon: pd.DataFrame, tfbs: pd.DataFrame) -> pd.DataFrame:

        tfbs_by_operon = tfbs.explode('tfbs_operon')
        operon = pd.merge(operon, tfbs_by_operon, how='left', left_on='operon_id', right_on='tfbs_operon')

        aggregation = {'tfbs_protrend_id': to_set_list, 'tfbs_operon': to_set_list}
        operon = self.group_by(df=operon, column='operon_id', aggregation=aggregation, default=flatten_set_list)

        operon = operon.rename(columns={'tfbs_protrend_id': 'tfbss'})
        operon = operon.drop(columns=['tfbs_operon'])

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
        operon = read_from_stack(stack=self.transform_stack, file='operon',
                                 default_columns=self.read_columns, reader=read_json_lines)
        non_empty_operon = operon['operon_id'] != ''
        operon = operon[non_empty_operon]

        operon = operon.explode('regulon')
        operon = self.drop_duplicates(df=operon, subset=['operon_id', 'regulon'], perfect_match=True, preserve_nan=True)

        aggregation = {'regulon': to_set_list}
        operon = self.group_by(df=operon, column='operon_id', aggregation=aggregation, default=flatten_set_list)

        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = self.select_columns(gene, 'protrend_id', 'locus_tag', 'name', 'locus_tag_old', 'strand', 'start', 'stop',
                                   'organism_protrend_id')

        gene = gene.dropna(subset=['protrend_id', 'locus_tag_old'])
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id',
                                    'locus_tag': 'gene_locus_tag',
                                    'name': 'gene_name',
                                    'locus_tag_old': 'gene_old_locus_tag',
                                    'strand': 'gene_strand',
                                    'start': 'gene_start',
                                    'stop': 'gene_stop'})

        tfbs = read_from_stack(stack=self.transform_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = self.select_columns(tfbs, 'protrend_id', 'operon')
        tfbs = tfbs.dropna(subset=['protrend_id'])
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id', 'operon': 'tfbs_operon'})

        # genes
        # tfbss
        # strand
        # start
        # stop
        operon = self._transform_operon_by_tfbs(operon=operon, tfbs=tfbs)
        operon = self._transform_operon_by_gene(operon=operon, gene=gene)
        operon = self._operon_coordinates(operon=operon)

        self._stack_transformed_nodes(operon)

        return operon


class OperonToGeneConnector(CollectfConnector):
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


class OperonToTFBSConnector(CollectfConnector):
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


class GeneToTFBSConnector(CollectfConnector):
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
