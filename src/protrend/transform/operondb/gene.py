from collections import defaultdict

import pandas as pd
from neo4j.exceptions import Neo4jError, DriverError

from protrend.io import read_txt, read
from protrend.model import Gene, Organism
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.operondb.base import OperonDBTransformer, OperonDBConnector
from protrend.transform.transformations import (drop_empty_string, merge_columns, create_input_value, group_by,
                                                select_columns)
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, to_int_str, flatten_set_list_nan, rstrip, lstrip, to_set_list


def _get_organism_from_gene(gene_node):
    organism_relationships = gene_node.organism.all()

    if organism_relationships:
        organism = organism_relationships[0]

        if organism.ncbi_taxonomy:
            return organism.protrend_id, organism.ncbi_taxonomy

        return organism.protrend_id, organism.name

    return None, None


class GeneTransformer(GeneMixIn, OperonDBTransformer,
                      source='operondb',
                      version='0.0.0',
                      node=Gene,
                      order=100,
                      register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop',
                       'organism', 'ncbi_taxonomy',
                       'operon_db_id'])

    def fetch_gene(self) -> pd.DataFrame:
        cols = ['locus_tag', 'ncbi_taxonomy', 'organism']
        try:
            genes = defaultdict(list)
            nodes = self.node.nodes.all()
            for node in nodes:
                val = getattr(node, 'locus_tag')
                genes['locus_tag'].append(val)

                organism_id, organism_taxonomy = _get_organism_from_gene(node)
                genes['organism'].append(organism_id)
                genes['ncbi_taxonomy'].append(organism_taxonomy)

            df = pd.DataFrame.from_dict(genes)

            if df.empty:
                return pd.DataFrame(columns=cols)

            return df

        except (Neo4jError, DriverError):
            return pd.DataFrame(columns=cols)

    @staticmethod
    def merge_operon_gene(operon: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        operon = operon.assign(locus_tag=operon['operon_genes'].copy())
        operon = operon.explode('locus_tag')
        operon = operon.assign(locus_tag=operon['locus_tag'].str.lower())

        gene = gene.assign(locus_tag=gene['locus_tag'].str.lower())

        operon_gene = pd.merge(operon, gene, on='locus_tag')
        operon_gene = operon_gene.drop(columns=['locus_tag'])
        operon_gene = operon_gene.rename(columns={'operon_genes': 'locus_tag'})

        aggregation = {'locus_tag': flatten_set_list_nan}
        operon_gene = group_by(operon_gene, column='operon_db_id', aggregation=aggregation)
        return operon_gene

    @staticmethod
    def transform_gene(operon_gene: pd.DataFrame) -> pd.DataFrame:
        genes = operon_gene.explode('locus_tag')
        genes = apply_processors(genes, locus_tag=[lstrip, rstrip], ncbi_taxonomy=to_int_str)

        aggregation = {'operon_db_id': to_set_list}
        genes = group_by(genes, column='locus_tag', aggregation=aggregation)
        genes = genes.dropna(subset=['locus_tag'])
        genes = drop_empty_string(genes, 'locus_tag')

        genes = create_input_value(genes, 'locus_tag')
        return genes

    def transform(self):
        conserved_columns = ['coid', 'org', 'name', 'op', 'definition', 'source', 'mbgd']
        known_columns = ['koid', 'org', 'name', 'op', 'definition', 'source']
        conserved = read(source=self.source, version=self.version, file='conserved_operon.txt', reader=read_txt,
                         default=pd.DataFrame(columns=conserved_columns))
        known = read(source=self.source, version=self.version, file='known_operon.txt', reader=read_txt,
                     default=pd.DataFrame(columns=known_columns))

        operon = self.transform_operon(conserved, known)
        operon = select_columns(operon, 'operon_db_id', 'operon_genes')

        gene = self.fetch_gene()

        operon_gene = self.merge_operon_gene(operon, gene)

        genes = self.transform_gene(operon_gene)
        annotated_genes = self.annotate_genes(genes)

        df = pd.merge(annotated_genes, genes, on='input_value', suffixes=('_annotation', '_operondb'))

        # merge loci
        df = merge_columns(df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_operondb')

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df


class GeneToOrganismConnector(OperonDBConnector,
                              source='operondb',
                              version='0.0.0',
                              from_node=Gene,
                              to_node=Organism,
                              register=True):

    def connect(self):
        df = self.create_connection(source='gene', target='gene',
                                    target_column='organism')
        self.stack_connections(df)
