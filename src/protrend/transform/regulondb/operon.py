from statistics import mode, StatisticsError

import numpy as np
import pandas as pd

from protrend.io import read_from_stack, read_json_frame, read_txt
from protrend.model.model import Operon, Gene, TFBS, Promoter
from protrend.utils.processors import (apply_processors, operon_name, flatten_set_list, to_list,
                                       to_set_list, operon_hash)
from protrend.transform.regulondb.base import RegulondbTransformer, RegulondbConnector
from protrend.transform.regulondb.gene import GeneTransformer
from protrend.transform.regulondb.promoter import PromoterTransformer
from protrend.transform.regulondb.tfbs import TFBSTransformer
from protrend.utils import SetList
from protrend.utils.miscellaneous import is_null


class OperonTransformer(RegulondbTransformer):
    default_node = Operon
    default_transform_stack = {'tu': 'transcription_unit.txt', 'tu_gene': 'tu_gene_link.txt',
                               'gene': 'integrated_gene.json', 'tfbs': 'integrated_tfbs.json',
                               'promoter': 'integrated_promoter.json'}
    default_order = 80
    columns = SetList(['name', 'promoters', 'genes', 'tfbss', 'strand', 'start', 'stop', 'operon_hash', 'protrend_id',
                       'transcription_unit_id', 'promoter_id', 'operon_id', 'gene_id',
                       'gene_name', 'gene_strand', 'gene_start', 'gene_stop'])

    tu_columns = SetList(['transcription_unit_id', 'promoter_id', 'transcription_unit_name', 'operon_id',
                          'key_id_org', 'transcription_unit_note', 'tu_internal_comment'])

    tu_gene_columns = SetList(['transcription_unit_id', 'gene_id'])

    def _transform_operon_by_gene(self, operon: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        operon = operon.explode(column='gene_id')

        operon = pd.merge(operon, gene, on='gene_id')

        # group by the operon_id
        aggregation = {'gene_id': to_set_list,
                       'gene_protrend_id': to_set_list,
                       'gene_name': to_set_list,
                       'gene_strand': to_set_list,
                       'gene_start': to_set_list,
                       'gene_stop': to_set_list}
        operon = self.group_by(df=operon, column='operon_id', aggregation=aggregation, default=flatten_set_list)

        operon = operon.rename(columns={'gene_protrend_id': 'genes'})

        operon['operon_hash'] = operon['genes']
        operon['name'] = operon['gene_name']

        operon = apply_processors(operon, operon_hash=[to_list, operon_hash], name=[to_list, operon_name])
        operon = operon.dropna(subset=['operon_hash'])
        operon = self.drop_duplicates(df=operon, subset=['operon_hash'], perfect_match=True, preserve_nan=True)

        return operon

    def _transform_operon_by_tfbs(self, operon: pd.DataFrame, tfbs: pd.DataFrame) -> pd.DataFrame:
        operon_tfbs = self._build_operon_tfbs()
        operon_tfbs = pd.merge(operon_tfbs, tfbs, on='site_id')

        operon = pd.merge(operon, operon_tfbs, how='left', on='operon_id')

        aggregation = {'tfbs_protrend_id': to_set_list, 'site_id': to_set_list}
        operon = self.group_by(df=operon, column='operon_id', aggregation=aggregation, default=flatten_set_list)

        operon = operon.rename(columns={'tfbs_protrend_id': 'tfbss'})

        return operon

    def _transform_operon_by_promoter(self, operon: pd.DataFrame, promoter: pd.DataFrame) -> pd.DataFrame:
        operon = operon.explode(column='promoter_id')
        operon = pd.merge(operon, promoter, how='left', on='promoter_id')

        aggregation = {'promoter_protrend_id': to_set_list, 'promoter_id': to_set_list}
        operon = self.group_by(df=operon, column='operon_id', aggregation=aggregation, default=flatten_set_list)

        operon = operon.rename(columns={'promoter_protrend_id': 'promoters'})

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
        tu = read_from_stack(stack=self.transform_stack, file='tu',
                             default_columns=self.tu_columns, reader=read_txt,
                             skiprows=34, names=self.tu_columns)

        tu_gene = read_from_stack(stack=self.transform_stack, file='tu_gene',
                                  default_columns=self.tu_gene_columns, reader=read_txt,
                                  skiprows=29, names=self.tu_gene_columns)

        tu = pd.merge(tu, tu_gene, on='transcription_unit_id')
        tu = self.select_columns(tu, 'transcription_unit_id', 'promoter_id', 'operon_id', 'gene_id')

        operon = self.group_by(tu, column='operon_id', aggregation={}, default=to_set_list)

        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = self.select_columns(gene, 'protrend_id', 'name', 'strand', 'start', 'stop', 'gene_id')
        gene = gene.dropna(subset=['protrend_id', 'gene_id'])
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id',
                                    'name': 'gene_name',
                                    'strand': 'gene_strand',
                                    'start': 'gene_start',
                                    'stop': 'gene_stop'})

        tfbs = read_from_stack(stack=self.transform_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = self.select_columns(tfbs, 'protrend_id', 'site_id')
        tfbs = tfbs.dropna(subset=['protrend_id', 'site_id'])
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id'})

        promoter = read_from_stack(stack=self.transform_stack, file='promoter',
                                   default_columns=PromoterTransformer.columns, reader=read_json_frame)
        promoter = self.select_columns(promoter, 'protrend_id', 'promoter_id')
        promoter = promoter.dropna(subset=['protrend_id', 'promoter_id'])
        promoter = promoter.rename(columns={'protrend_id': 'promoter_protrend_id'})

        # genes
        # promoters
        # tfbss
        # strand
        # start
        # stop
        operon = self._transform_operon_by_promoter(operon=operon, promoter=promoter)
        operon = self._transform_operon_by_tfbs(operon=operon, tfbs=tfbs)
        operon = self._transform_operon_by_gene(operon=operon, gene=gene)
        operon = self._operon_coordinates(operon=operon)

        self._stack_transformed_nodes(operon)

        return operon


class OperonToGeneConnector(RegulondbConnector):
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


class OperonToPromoterConnector(RegulondbConnector):
    default_from_node = Operon
    default_to_node = Promoter
    default_connect_stack = {'operon': 'integrated_operon.json'}

    def connect(self):
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)

        operon = apply_processors(operon, promoters=to_list)
        operon = operon.explode('promoters')

        operon = operon.dropna(subset=['protrend_id', 'promoters'])
        operon = operon.drop_duplicates(subset=['protrend_id', 'promoters'])

        from_identifiers = operon['protrend_id'].tolist()
        to_identifiers = operon['promoters'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class OperonToTFBSConnector(RegulondbConnector):
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


class GeneToTFBSConnector(RegulondbConnector):
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


class GeneToPromoterConnector(RegulondbConnector):
    default_from_node = Gene
    default_to_node = Promoter
    default_connect_stack = {'operon': 'integrated_operon.json'}

    def connect(self):
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)

        operon = apply_processors(operon, genes=to_list, promoters=to_list)
        operon = operon.explode('promoters')
        operon = operon.explode('genes')

        operon = operon.dropna(subset=['genes', 'promoters'])
        operon = operon.drop_duplicates(subset=['genes', 'promoters'])

        from_identifiers = operon['genes'].tolist()
        to_identifiers = operon['promoters'].tolist()
        kwargs = dict(operon=operon['protrend_id'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_json(df)
