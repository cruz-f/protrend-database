import pandas as pd

from protrend.io.csv import read_csv
from protrend.io.json import read_json_lines
from protrend.transform.connector import DefaultConnector
from protrend.transform.processors import apply_processors, str_join, operon_name, genes_to_hash, operon_strand, \
    operon_left_position, operon_right_position
from protrend.transform.regprecise.gene import GeneTransformer
from protrend.transform.regprecise.settings import OperonSettings, OperonToSource
from protrend.transform.regprecise.source import SourceTransformer
from protrend.transform.regprecise.tfbs import TFBSTransformer
from protrend.transform.transformer import DefaultTransformer
from protrend.utils.graph import build_graph, find_connected_nodes
from protrend.utils.miscellaneous import flatten_list


class OperonTransformer(DefaultTransformer):
    default_settings = OperonSettings
    columns = {'protrend_id',
               'operon_id', 'name', 'url', 'regulon', 'tfbs',
               'tfbss',
               'operon_id_old', 'operon_id_new', 'locus_tag', 'locus_tag_regprecise',
               'genes',
               'first_gene_position_left',
               'last_gene_position_right', }

    def _read_operon(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('operon')

        if file_path:
            df = read_json_lines(file_path)

        else:
            df = pd.DataFrame(columns=['operon_id', 'name', 'url', 'regulon', 'tfbs', 'gene'])

        return df

    def _read_gene(self) -> pd.DataFrame:

        file_path = self._transform_stack.get('gene')

        if file_path:
            df = read_csv(file_path)

        else:
            df = pd.DataFrame(columns=[GeneTransformer.columns])

        df = df.rename(columns={'protrend_id': 'gene_protrend_id'})

        df = df.dropna(subset=['gene_protrend_id'])

        df = df[['gene_protrend_id', 'locus_tag', 'name', 'locus_tag_regprecise',
                 'strand', 'position_left', 'position_right']]

        return df

    def _read_tfbs(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('tfbs')

        if file_path:
            df = read_csv(file_path)

        else:
            df = pd.DataFrame(columns=[TFBSTransformer.columns])

        df = df.dropna(subset=['protrend_id'])

        df = df.rename(columns={'protrend_id': 'tfbs_protrend_id'})

        df = df[['tfbs_protrend_id', 'operon']]

        return df

    @staticmethod
    def _operon_by_gene(operon: pd.DataFrame) -> pd.DataFrame:

        # group duplicates
        operon_by_gene = operon.explode('gene')

        agg_funcs = {'operon_id': set,
                     'url': flatten_list,
                     'regulon': flatten_list,
                     'tfbs': flatten_list,
                     'tfbss': flatten_list}

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
                     'tfbs': flatten_list,
                     'tfbss': flatten_list,
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
                     'tfbs': flatten_list,
                     'tfbss': flatten_list,
                     'name': set,
                     'gene_protrend_id': set,
                     'locus_tag': set,
                     'locus_tag_regprecise': set}

        operon = operon_by_gene.groupby(['operon_id_new']).aggregate(agg_funcs)
        operon = operon.reset_index()
        operon = operon.rename(columns={'gene_protrend_id': 'genes'})

        operon['operon_id'] = operon['genes']

        apply_processors(str_join,
                         df=operon,
                         col='operon_id')

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
                     'tfbs': flatten_list,
                     'gene': flatten_list,
                     'tfbs_protrend_id': set}

        operon = operon.groupby(['operon_id']).aggregate(agg_funcs)
        operon = operon.reset_index()

        operon = operon.rename(columns={'tfbs_protrend_id': 'tfbss'})

        return operon

    @staticmethod
    def _operon_coordinates(operon: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:

        gene = gene.set_index(gene['gene_protrend_id'])

        strands = []
        positions_left = []
        positions_right = []

        for genes in operon['genes']:

            op_strand = None

            for op_gene in genes:
                op_gene_strand = gene.loc[op_gene, 'strand']
                op_strand = operon_strand(previous_strand=op_strand,
                                          current_strand=op_gene_strand)

            operon_left = None
            operon_right = None

            for op_gene in genes:
                op_gene_left = gene.loc[op_gene, 'position_left']
                op_gene_right = gene.loc[op_gene, 'position_right']

                operon_left = operon_left_position(strand=op_strand,
                                                   previous_left=operon_left,
                                                   current_left=op_gene_left)

                operon_right = operon_right_position(strand=op_strand,
                                                     previous_right=operon_right,
                                                     current_right=op_gene_right)

            strands.append(operon_strand)
            positions_left.append(positions_left)
            positions_right.append(positions_right)

        operon['first_gene_position_left'] = positions_left
        operon['last_gene_position_right'] = positions_right

        return operon

    def transform(self):
        operon = self._read_operon()
        operon = self.drop_duplicates(df=operon, subset=['operon_id'], perfect_match=True, preserve_nan=True)
        operon = operon.drop(columns=['name'], axis=1)

        gene = self._read_gene()
        tfbs = self._read_tfbs()

        # genes
        # tfbss
        # strand
        # first_gene_position_left
        # last_gene_position_right
        operon = self._transform_operon_by_tfbs(operon=operon, tfbs=tfbs)
        operon = self._transform_operon_by_gene(operon=operon, gene=gene)

        df = self._operon_coordinates(operon=operon, gene=gene)

        if df.empty:
            df = self.make_empty_frame()

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df

    def integrate(self, df: pd.DataFrame) -> pd.DataFrame:

        # ensure uniqueness
        df = self.drop_duplicates(df=df,
                                  subset=self.node_factors,
                                  perfect_match=False,
                                  preserve_nan=True)

        df['genes_id'] = df['genes'].map(genes_to_hash)

        # take a db snapshot for the current node
        snapshot = self.node_view()
        snapshot['genes_id'] = snapshot['genes'].map(genes_to_hash)

        # find matching nodes according to several node factors/properties
        nodes_mask = self.find_nodes(nodes=df, snapshot=snapshot, node_factors=('genes_id',))

        # nodes to be updated
        update_nodes = df[nodes_mask]

        # find/set protrend identifiers for update nodes
        update_size, _ = update_nodes.shape
        ids_mask = self.find_snapshot(nodes=update_nodes, snapshot=snapshot, node_factors=('genes_id',))
        update_nodes[self.node.identifying_property] = snapshot.loc[ids_mask, self.node.identifying_property]
        update_nodes['load'] = ['update'] * update_size
        update_nodes['what'] = ['nodes'] * update_size

        # nodes to be created
        create_nodes = df[~nodes_mask]

        # create/set new protrend identifiers
        create_size, _ = create_nodes.shape
        create_identifiers = self.protrend_identifiers_batch(create_size)
        create_nodes[self.node.identifying_property] = create_identifiers
        create_nodes['load'] = ['create'] * create_size
        update_nodes['what'] = ['nodes'] * update_size

        # concat both dataframes
        df = pd.concat([create_nodes, update_nodes], axis=0)
        df.drop(columns=['genes_id'], axis=1)
        df_name = f'integrated_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        self._stack_integrated_nodes(df)

        return df


class OperonToSourceConnector(DefaultConnector):
    default_settings = OperonToSource

    def _read_operon(self) -> pd.DataFrame:
        file_path = self._connect_stack.get('operon')

        if file_path:
            df = read_csv(file_path)

        else:
            df = pd.DataFrame(columns=OperonTransformer.columns)

        return df

    def _read_source(self) -> pd.DataFrame:
        file_path = self._connect_stack.get('source')

        if file_path:
            df = read_csv(file_path)

        else:
            df = pd.DataFrame(columns=SourceTransformer.columns)

        return df

    def connect(self):

        operon = self._read_operon()
        operon = operon.explode('regulon')
        source = self._read_source()

        from_identifiers = operon['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=operon['url'].tolist(),
                      external_identifier=operon['regulon'].tolist(),
                      key=['regulon_id'] * size)

        df = self.make_connection(size=size,
                                  from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_csv(df)
