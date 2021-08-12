import pandas as pd

from protrend.io.utils import read_from_stack
from protrend.transform.connector import DefaultConnector
from protrend.transform.processors import apply_processors, str_join, operon_name, genes_to_hash, operon_strand, \
    operon_left_position, operon_right_position
from protrend.transform.regprecise.gene import GeneTransformer
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.regprecise.settings import OperonSettings, OperonToSource, OperonToOrganism, OperonToRegulator, \
    OperonToGene, OperonToTFBS
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
    read_columns = {'operon_id', 'name', 'url', 'regulon', 'tfbs', 'gene'}

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
        operon = operon.merge(tfbs_by_operon, left_on='operon_id', right_on='operon')

        agg_funcs = {'url': set,
                     'regulon': set,
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
        operon = read_from_stack(tl=self, file='operon', json=True, default_columns=self.read_columns)
        operon = operon.drop(columns=['name'], axis=1)
        operon = operon.explode('regulon')
        operon = self.drop_duplicates(df=operon, subset=['operon_id', 'regulon'],
                                      perfect_match=True, preserve_nan=True)

        gene = read_from_stack(tl=self, file='gene', json=False, default_columns=GeneTransformer.columns)
        gene = gene[['protrend_id', 'locus_tag', 'name', 'locus_tag_regprecise',
                     'strand', 'position_left', 'position_right']]
        gene = gene.dropna(subset=['protrend_id'])
        gene = gene.dropna(subset=['locus_tag_regprecise'])
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id'})

        tfbs = read_from_stack(tl=self, file='tfbs', json=False, default_columns=TFBSTransformer.columns)
        tfbs = tfbs[['protrend_id', 'operon']]
        tfbs = tfbs.dropna(subset=['protrend_id'])
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id'})

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

        df['genes_id'] = df['genes'].map(genes_to_hash)

        # ensure uniqueness
        df = self.drop_duplicates(df=df,
                                  subset=('genes_id', ),
                                  perfect_match=True,
                                  preserve_nan=True)

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

    def connect(self):
        operon = read_from_stack(tl=self, file='operon', json=False, default_columns=OperonTransformer.columns)
        operon = operon.explode('regulon')
        source = read_from_stack(tl=self, file='source', json=False, default_columns=SourceTransformer.columns)

        from_identifiers = operon['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=operon['url'].tolist(),
                      external_identifier=operon['regulon'].tolist(),
                      key=['regulon_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_csv(df)


class OperonToOrganismConnector(DefaultConnector):
    default_settings = OperonToOrganism

    def connect(self):
        operon = read_from_stack(tl=self, file='operon', json=False, default_columns=OperonTransformer.columns)
        apply_processors(list, df=operon, col='regulon')
        operon = operon.explode('regulon')
        regulator = read_from_stack(tl=self, file='regulator', json=False, default_columns=RegulatorTransformer.columns)

        merged = pd.merge(operon, regulator, left_on='regulon', right_on='regulon_id',
                          suffixes=('_operon', '_regulator'))
        merged = merged.dropna(subset=['protrend_id_operon'])
        merged = merged.dropna(subset=['organism_protrend_id'])
        merged = merged.drop_duplicates(subset=['protrend_id_operon', 'organism_protrend_id'])

        from_identifiers = merged['protrend_id_operon'].tolist()
        to_identifiers = merged['organism_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)


class OperonToRegulatorConnector(DefaultConnector):
    default_settings = OperonToRegulator

    def connect(self):
        operon = read_from_stack(tl=self, file='operon', json=False, default_columns=OperonTransformer.columns)
        apply_processors(list, df=operon, col='regulon')
        operon = operon.explode('regulon')
        regulator = read_from_stack(tl=self, file='regulator', json=False, default_columns=RegulatorTransformer.columns)

        merged = pd.merge(operon, regulator, left_on='regulon', right_on='regulon_id',
                          suffixes=('_operon', '_regulator'))
        merged = merged.dropna(subset=['protrend_id_operon'])
        merged = merged.dropna(subset=['protrend_id_regulator'])
        merged = merged.drop_duplicates(subset=['protrend_id_operon', 'protrend_id_regulator'])

        from_identifiers = merged['protrend_id_operon'].tolist()
        to_identifiers = merged['protrend_id_regulator'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)


class OperonToGeneConnector(DefaultConnector):
    default_settings = OperonToGene

    def connect(self):
        operon = read_from_stack(tl=self, file='operon', json=False, default_columns=OperonTransformer.columns)
        apply_processors(list, df=operon, col='genes')
        operon = operon.explode('genes')
        operon = operon.dropna(subset=['protrend_id'])
        operon = operon.dropna(subset=['genes'])
        operon = operon.drop_duplicates(subset=['protrend_id', 'genes'])

        from_identifiers = operon['protrend_id'].tolist()
        to_identifiers = operon['genes'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)


class OperonToTFBSConnector(DefaultConnector):
    default_settings = OperonToTFBS

    def connect(self):
        operon = read_from_stack(tl=self, file='operon', json=False, default_columns=OperonTransformer.columns)
        apply_processors(list, df=operon, col='tfbss')
        operon = operon.explode('tfbss')
        operon = operon.dropna(subset=['protrend_id'])
        operon = operon.dropna(subset=['tfbss'])
        operon = operon.drop_duplicates(subset=['protrend_id', 'tfbss'])

        from_identifiers = operon['protrend_id'].tolist()
        to_identifiers = operon['tfbss'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)
