import pandas as pd

from protrend.io import read_json_frame, read_json_lines, read_from_stack
from protrend.model import TFBS
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.transform.dbtbs.gene import GeneTransformer
from protrend.utils.processors import (apply_processors, to_list, to_str, operon_hash, site_hash, to_set_list,
                                       take_first, upper_case, flatten_set_list)
from protrend.utils import SetList, is_null


class TFBSTransformer(DBTBSTransformer,
                      source='dbtbs',
                      version='0.0.3',
                      node=TFBS,
                      order=90,
                      register=True):
    default_transform_stack = {'tfbs': 'TFBS.json',
                               'gene': 'integrated_gene.json'}
    columns = SetList(['sequence', 'strand', 'start', 'stop', 'length', 'site_hash', 'protrend_id',
                       'identifier', 'url', 'regulation', 'pubmed', 'tf', 'operon', 'gene',
                       'location', 'absolute_position'])
    read_columns = SetList(['identifier', 'url', 'regulation', 'pubmed', 'tf', 'operon', 'gene',
                            'location', 'absolute_position', 'sequence'])

    def _transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        # filter by nan
        tfbs = tfbs.dropna(subset=['identifier', 'sequence'])

        # filter by id
        tfbs = self.drop_duplicates(df=tfbs, subset=['identifier'], perfect_match=True)

        # processing
        tfbs = tfbs.explode(column='sequence')
        tfbs = apply_processors(tfbs, sequence=upper_case, location=take_first, absolute_position=take_first)

        # filter by coordinates
        tfbs = self.drop_duplicates(df=tfbs, subset=['sequence', 'location', 'absolute_position'], perfect_match=True)

        return tfbs

    @staticmethod
    def _tfbs_coordinates(tfbs: pd.DataFrame) -> pd.DataFrame:

        def sequence_len(item):
            if is_null(item):
                return None

            if isinstance(item, str):
                return len(item)

            return None

        def sequence_start(item):
            if is_null(item):
                return None

            if isinstance(item, str):
                try:
                    start, _ = item.split('..')

                    start = int(start)
                    return start

                except ValueError:
                    return None

            return None

        def sequence_stop(item):
            if is_null(item):
                return None

            if isinstance(item, str):
                try:
                    _, stop = item.split('..')

                    stop = int(stop)
                    return stop

                except ValueError:
                    return None

            return None

        def sequence_strand(item):
            if is_null(item):
                return None

            if isinstance(item, str):
                try:
                    start, stop = item.split('..')

                    start = int(start)
                    stop = int(stop)
                    diff = stop - start

                    if diff < 0:
                        return 'reverse'

                    return 'forward'

                except ValueError:
                    return None

            return None

        tfbs['length'] = tfbs['sequence'].map(sequence_len, na_action='ignore')
        tfbs['start'] = tfbs['absolute_position'].map(sequence_start, na_action='ignore')
        tfbs['stop'] = tfbs['absolute_position'].map(sequence_stop, na_action='ignore')
        tfbs['strand'] = tfbs['absolute_position'].map(sequence_strand, na_action='ignore')

        return tfbs

    def transform(self):
        tfbs = read_from_stack(stack=self.transform_stack, key='tfbs', columns=self.read_columns,
                               reader=read_json_lines)

        tfbs = self._transform_tfbs(tfbs)
        tfbs = self._tfbs_coordinates(tfbs)

        gene = read_from_stack(stack=self.transform_stack, key='gene', columns=GeneTransformer.columns,
                               reader=read_json_frame)
        gene = self.select_columns(gene, 'protrend_id', 'name_dbtbs')
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id', 'name_dbtbs': 'gene_name_dbtbs'})

        tfbs = apply_processors(tfbs, gene=flatten_set_list)
        tfbs = tfbs.explode(column='gene')

        gene_tfbs = pd.merge(tfbs, gene, left_on='gene', right_on='gene_name_dbtbs')
        aggregation = {'location': take_first,
                       'absolute_position': take_first,
                       'sequence': take_first,
                       'length': take_first,
                       'start': take_first,
                       'stop': take_first,
                       'strand': take_first,
                       'gene': to_set_list,
                       'gene_name_dbtbs': to_set_list,
                       'gene_protrend_id': to_set_list}

        gene_tfbs = self.group_by(df=gene_tfbs, column='identifier', aggregation=aggregation, default=flatten_set_list)

        # filter by site hash: length + strand + start + genes
        df = apply_processors(gene_tfbs,
                              sequence=[to_str, to_list],
                              length=[to_str, to_list],
                              strand=[to_str, to_list],
                              start=[to_str, to_list],
                              gene_protrend_id=[to_list, operon_hash, to_list])

        gene_tfbs['site_hash'] = df['sequence'] + df['length'] + df['strand'] + df['start'] + df['gene_protrend_id']
        gene_tfbs = apply_processors(gene_tfbs, site_hash=site_hash)

        gene_tfbs = self.drop_duplicates(df=gene_tfbs, subset=['site_hash'], perfect_match=True)
        gene_tfbs = gene_tfbs.dropna(subset=['site_hash'])

        self.stack_transformed_nodes(gene_tfbs)

        return gene_tfbs
