import pandas as pd

from protrend.io.json import read_json_lines, read_json_frame
from protrend.io.utils import read_from_stack
from protrend.model.model import TFBS
from protrend.transform.collectf.base import CollectfTransformer
from protrend.transform.collectf.gene import GeneTransformer
from protrend.transform.processors import (apply_processors, to_list, flatten_set_list,
                                           take_last, genes_to_hash, to_str, to_set_list)
from protrend.utils.miscellaneous import is_null


class TFBSTransformer(CollectfTransformer):
    default_node = TFBS
    default_node_factors = ('site_hash', )
    default_transform_stack = {'tfbs': 'TFBS.json', 'gene': 'integrated_gene.json'}
    default_order = 70
    columns = {'protrend_id',
               'sequence', 'strand', 'start', 'stop', 'length', 'site_hash',
               'tfbs_id', 'mode', 'pubmed', 'organism', 'regulon', 'operon', 'gene', 'experimental_evidence'}
    read_columns = {'tfbs_id', 'site_start', 'site_end', 'site_strand', 'mode', 'sequence', 'pubmed', 'organism',
                    'regulon', 'operon', 'gene', 'experimental_evidence'}

    def _transform_tfbs(self, tfbs: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:

        # filter by sequence
        tfbs = tfbs.dropna(subset=['sequence'])

        # filter by gene
        tfbs = apply_processors(tfbs, gene=to_list)
        tfbs = tfbs.explode(column='gene')

        tfbs = pd.merge(tfbs, gene, left_on='gene', right_on='gene_old_locus_tag')

        aggr = {'pubmed': flatten_set_list, 'regulon': flatten_set_list, 'operon': flatten_set_list,
                'experimental_evidence': flatten_set_list, 'gene': to_set_list}
        tfbs = self.group_by(df=tfbs, column='tfbs_id', aggregation=aggr, default=take_last)

        # filter by regulon, sequence and position
        tfbs = apply_processors(tfbs, regulon=to_list)
        tfbs = tfbs.explode(column='regulon')
        tfbs = self.drop_duplicates(df=tfbs, subset=['tfbs_id', 'sequence', 'regulon'],
                                    perfect_match=True, preserve_nan=True)
        tfbs = tfbs.reset_index(drop=True)

        return tfbs

    @staticmethod
    def _tfbs_coordinates(tfbs: pd.DataFrame) -> pd.DataFrame:

        def strand_translation(item):
            if is_null(item):
                return None

            if item == '1':
                return 'forward'

            elif item == '-1':
                return 'reverse'

            return None

        def sequence_len(item):
            if is_null(item):
                return None

            if isinstance(item, str):
                return len(item)

            return None

        tfbs['length'] = tfbs['sequence'].map(sequence_len, na_action='ignore')
        tfbs['strand'] = tfbs['strand'].map(strand_translation, na_action='ignore')

        return tfbs

    def transform(self):
        tfbs = read_from_stack(stack=self.transform_stack, file='tfbs',
                               default_columns=self.read_columns, reader=read_json_lines)
        tfbs = tfbs.rename(columns={'site_start': 'start', 'site_end': 'stop', 'site_strand': 'strand'})

        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = self.select_columns(gene, 'protrend_id', 'locus_tag_old')
        gene = gene.rename(columns={'locus_tag_old': 'gene_old_locus_tag', 'protrend_id': 'gene_protrend_id'})
        gene = gene.dropna(subset=['gene_old_locus_tag', 'gene_protrend_id'])
        gene = self.drop_duplicates(df=gene, subset=['gene_old_locus_tag', 'gene_protrend_id'],
                                    perfect_match=True, preserve_nan=False)

        df = self._transform_tfbs(tfbs=tfbs, gene=gene)

        df = self._tfbs_coordinates(df)

        # filter by site hash: length + strand + start + genes
        str_sequence = df['sequence'].map(to_str, na_action='ignore').fillna('')
        str_len = df['length'].map(to_str, na_action='ignore').fillna('')
        str_strand = df['strand'].map(to_str, na_action='ignore').fillna('')
        str_start = df['start'].map(to_str, na_action='ignore').fillna('')
        str_genes = df['gene_protrend_id'].map(genes_to_hash, na_action='ignore').fillna('')
        df['site_hash'] = str_sequence + str_len + str_strand + str_start + str_genes

        df = self.drop_duplicates(df=df, subset=['site_hash'], perfect_match=True, preserve_nan=True)

        self._stack_transformed_nodes(df)

        return df
