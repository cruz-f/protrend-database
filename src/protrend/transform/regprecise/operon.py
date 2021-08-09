from collections import defaultdict
from typing import List, Union

import pandas as pd

from protrend.io.csv import read_csv
from protrend.io.json import read_json_lines
from protrend.model.model import Source, Organism
from protrend.transform.annotation.gene import annotate_genes
from protrend.transform.dto import GeneDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, nan_to_str, take_first, null_to_nan, \
    str_join, operon_name
from protrend.transform.regprecise.settings import GeneSettings, OperonSettings
from protrend.transform.transformer import Transformer
from protrend.utils.graph import build_graph, find_connected_nodes
from protrend.utils.miscellaneous import take_last, flatten_list


class OperonTransformer(Transformer):

    def __init__(self, settings: OperonSettings = None):

        if not settings:
            settings = OperonSettings()

        super().__init__(settings)

    def _read_operon(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('operon')

        if not file_path:
            return pd.DataFrame(columns=['operon_id', 'name', 'url', 'regulon', 'tfbs', 'gene'])

        return read_json_lines(file_path)

    def _read_gene(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('gene')

        if not file_path:
            return pd.DataFrame(columns=['gene_protrend_id', 'locus_tag', 'name', 'locus_tag_regprecise'])

        df = read_csv(file_path)

        df = df.dropna(subset=['protrend_id'])

        df = df.rename(columns={'protrend_id': 'gene_protrend_id'})

        df = df[['gene_protrend_id', 'locus_tag', 'name', 'locus_tag_regprecise']]

        return df

    def _read_tfbs(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('tfbs')

        if not file_path:
            return pd.DataFrame(columns=['tfbs_protrend_id', 'tfbs_id'])

        df = read_csv(file_path)

        df = df.dropna(subset=['tfbs_protrend_id'])

        df = df.rename(columns={'protrend_id': 'tfbs_protrend_id'})

        df = df[['tfbs_protrend_id', 'tfbs_id', 'operon']]

        return df

    @staticmethod
    def _operon_by_gene(operon: pd.DataFrame) -> pd.DataFrame:

        # group duplicates
        operon_by_gene = operon.explode('gene')
        agg_funcs = {'operon_id': set,
                     'url': set,
                     'regulon': flatten_list,
                     'tfbss': flatten_list,
                     'tfbs_old': flatten_list}
        operon_by_gene = operon_by_gene.groupby(['gene']).aggregate(agg_funcs)
        operon_by_gene = operon_by_gene.reset_index()

        return operon_by_gene

    @staticmethod
    def _normalize_operon(genes) -> pd.DataFrame:

        graph = build_graph(genes)
        operons = find_connected_nodes(graph)

        operon = pd.DataFrame({'gene': operons})

        apply_processors(list,
                         df=operon,
                         col='gene')

        operon['operon_id'] = operon['gene']

        apply_processors(str_join,
                         df=operon,
                         col='operon_id')

        return operon

    def _transform_operon_by_gene(self, operon: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:

        genes = operon['gene'].tolist()
        normalized_operon = self._normalize_operon(genes)
        normalized_operon_by_gene = normalized_operon.explode('gene')

        operon_by_gene = self._operon_by_gene(operon)

        operon = pd.merge(normalized_operon_by_gene, operon_by_gene, on='gene', suffixes=("_new", "_old"))

        # group by the new operons
        agg_funcs = {'operon_id_old': flatten_list,
                     'url': flatten_list,
                     'regulon': flatten_list,
                     'tfbss': flatten_list,
                     'tfbs_old': flatten_list,
                     'gene': set}

        operon = operon.groupby(['operon_id_new']).aggregate(agg_funcs)
        operon = operon.reset_index()

        apply_processors(list,
                         df=operon,
                         col='gene')

        operon_by_gene = operon.explode('gene')
        operon_by_gene = operon_by_gene.merge(gene, left_on='gene', right_on='locus_tag_regprecise')

        agg_funcs = {'operon_id_old': flatten_list,
                     'url': flatten_list,
                     'regulon': flatten_list,
                     'tfbss': flatten_list,
                     'tfbs_old': flatten_list,
                     'gene_protrend_id': set,
                     'name': set}

        operon = operon_by_gene.groupby(['operon_id_new']).aggregate(agg_funcs)
        operon = operon.reset_index()
        operon = operon.rename(columns={'gene_protrend_id': 'genes'})

        operon['operon_id_new'] = operon['genes']

        apply_processors(str_join,
                         df=operon,
                         col='operon_id_new')

        apply_processors(operon_name,
                         df=operon,
                         col='name')

        return operon

    @staticmethod
    def _transform_operon_by_tfbs(operon: pd.DataFrame, tfbs: pd.DataFrame) -> pd.DataFrame:

        tfbs_by_operon = tfbs.explode('operon')
        operon = operon.merge(tfbs_by_operon, how='inner', left_on='operon_id', right_on='operon')

        agg_funcs = {'url': set,
                     'regulon': flatten_list,
                     'tfbs_id': set,
                     'tfbs_protrend_id': set,
                     'gene': flatten_list}

        operon = operon.groupby(['operon_id']).aggregate(agg_funcs)

        operon = operon.rename(columns={'tfbs_id': 'tfbs_old',
                                        'tfbs_protrend_id': 'tfbss'})

        return operon

    def _add_operon_coordinates(self, operon: pd.DataFrame) -> pd.DataFrame:
        file_path = self._transform_stack.get('gene')

        if not file_path:
            df = pd.DataFrame(columns=['gene_protrend_id', 'strand', 'position_left', 'position_right'])

        else:
            df = read_csv(file_path)

            df = df.dropna(subset=['protrend_id'])

            df = df.rename(columns={'protrend_id': 'gene_protrend_id'})

        df = df.set_index(df['gene_protrend_id'])

        strands = []
        positions_left = []
        positions_right = []

        for genes in operon['genes']:

            operon_strand = None
            operon_left = None
            operon_right = None

            for gene in genes:

                gene_strand = df.loc[gene, 'strand']
                gene_left = df.loc[gene, 'position_left']
                gene_right = df.loc[gene, 'position_right']

                if gene_strand and operon_strand is None:
                    operon_strand = str(gene_strand)

                if gene_left and operon_left is None:
                    operon_left = int(gene_left)

                elif gene_left and gene_left < operon_left:
                    operon_left = int(gene_left)

                if gene_right and operon_right is None:
                    operon_right = int(gene_right)

                elif gene_right and gene_right > operon_right:
                    operon_right = int(gene_right)

            if operon_strand is None:
                operon_strand = 'forward'

            strands.append(operon_strand)
            positions_left.append(positions_left)
            positions_right.append(positions_right)

        operon['first_gene_position_left'] = positions_left
        operon['last_gene_position_right'] = positions_right

        return operon

    def transform(self):
        operon = self._read_operon()
        operon = self.drop_duplicates(df=operon, subset=['operon_id'], perfect_match=True, preserve_nan=True)
        operon = operon.drop(['name'], axis=1)

        gene = self._read_gene()
        tfbs = self._read_tfbs()

        # genes
        # tfbss
        # strand
        # first_gene_position_left
        # last_gene_position_right
        operon = self._transform_operon_by_tfbs(operon=operon, tfbs=tfbs)
        operon = self._transform_operon_by_gene(operon=operon, gene=gene)
        operon = self._add_operon_coordinates(operon=operon)

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, operon)

        return operon

    def _connect_to_source(self) -> pd.DataFrame:

        from_path = self._connect_stack.get('from')
        to_path = self._connect_stack.get('to_source')

        if not from_path:
            return pd.DataFrame()

        if not to_path:
            return pd.DataFrame()

        from_df = read_csv(from_path)

        from_identifiers = []
        kwargs = defaultdict(list)

        for _, row in from_df.iterrows():

            from_id = row['protrend_id']
            urls = row['url']
            regulons = row['regulon']

            for url, regulon in zip(urls, regulons):

                from_identifiers.append(from_id)
                kwargs['url'].append(url)
                kwargs['external_identifier'].append(regulon)
                kwargs['key'].append('regulon_id')

        size = len(from_identifiers)

        to_df = read_csv(to_path)
        to_df = to_df.query('name == regprecise')
        regprecise_id = to_df['protrend_id'].iloc[0]
        to_identifiers = [regprecise_id] * size

        return self.make_connection(size=size,
                                    from_node=self.node,
                                    to_node=Source,
                                    from_identifiers=from_identifiers,
                                    to_identifiers=to_identifiers,
                                    kwargs=kwargs)

    def connect(self):
        connection = self._connect_to_source()
        df_name = f'connected_{self.node.node_name()}_{Source.node_name()}'
        self.stack_csv(df_name, connection)

        #
        # connection = self._connect_to_organism()
        # df_name = f'connected_{self.node.node_name()}_{Organism.node_name()}'
        # self.stack_csv(df_name, connection)
