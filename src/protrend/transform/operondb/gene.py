from collections import defaultdict
from typing import Tuple

import pandas as pd
from neo4j.exceptions import Neo4jError

from protrend.io import read_from_stack, read_txt
from protrend.model import Gene, Organism
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.operondb.base import OperonDBTransformer, OperonDBConnector
from protrend.transform.transformations import drop_empty_string, drop_duplicates, merge_columns, create_input_value
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, to_list_nan


class GeneTransformer(GeneMixIn, OperonDBTransformer,
                      source='operondb',
                      version='0.0.0',
                      node=Gene,
                      order=100,
                      register=True):
    default_transform_stack = {'conserved': 'conserved_operon.txt', 'known': 'known_operon.txt'}

    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop',
                       'organism', 'ncbi_taxonomy',
                       'operon_db_id', 'operon_name', 'operon_function', 'pubmed'])
    conserved_columns = SetList(['coid', 'org', 'name', 'op', 'definition', 'source', 'mbgd'])
    known_columns = SetList(['koid', 'org', 'name', 'op', 'definition', 'source'])

    @staticmethod
    def transform_operon(conserved: pd.DataFrame, known: pd.DataFrame) -> pd.DataFrame:
        conserved = conserved.dropna(subset=['coid', 'op'])
        conserved = drop_empty_string(conserved, 'coid', 'op')
        conserved = drop_duplicates(conserved, subset=['coid', 'op'], perfect_match=True)

        known = known.dropna(subset=['koid', 'op'])
        known = drop_empty_string(known, 'koid', 'op')
        known = drop_duplicates(known, subset=['koid', 'op'], perfect_match=True)

        conserved = conserved.assign(operon_db_id=conserved['coid'].copy(),
                                     locus_tag=conserved['op'].str.split(','))
        known = known.assign(operon_db_id=known['koid'],
                             locus_tag=known['op'].str.split(','))

        operon = pd.concat([conserved, known])
        operon = operon.reset_index(drop=True)

        operon = operon.drop(columns=['koid', 'coid', 'op', 'mbgd'])
        operon = operon.rename(columns={'name': 'operon_name', 'definition': 'operon_function',
                                        'source': 'pubmed', 'org': 'ncbi_taxonomy'})
        return operon

    def fetch_gene(self) -> pd.DataFrame:

        def get_organism(gene_node):
            organism_relationships = gene_node.organism.all()

            if organism_relationships:
                organism = organism_relationships[0]

                if organism.ncbi_taxonomy:
                    return organism.protrend_id, organism.ncbi_taxonomy

                return organism.protrend_id, organism.name

            return None, None

        try:
            genes = defaultdict(list)
            nodes = self.node.nodes.all()
            for node in nodes:
                node_properties = node.properties

                for key in self.node.node_keys():
                    val = node_properties.get(key, None)
                    genes[key].append(val)

                organism_id, organism_taxonomy = get_organism(node)
                genes['organism'].append(organism_id)
                genes['ncbi_taxonomy'].append(organism_taxonomy)

            genes = pd.DataFrame.from_dict(genes)

        except Neo4jError:
            genes = pd.DataFrame(columns=list(self.node.node_keys()))
        return genes

    @staticmethod
    def merge_operon_gene(operon: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        operon = apply_processors(operon, locus_tag=to_list_nan)
        operon = operon.explode('locus_tag')

        operon = operon.assign(locus_tag_lower=operon['locus_tag'].str.lower())
        gene = gene.assign(locus_tag_lower=gene['locus_tag'].str.lower())

        operon = pd.merge(operon, gene, how='left', on='locus_tag_lower', suffixes=('_operon', '_gene'))

        operon = operon.drop(columns=['locus_tag_lower'])

        operon = merge_columns(operon, column='locus_tag', right='locus_tag_gene',
                               left='locus_tag_operon')
        operon = merge_columns(operon, column='ncbi_taxonomy', left='ncbi_taxonomy_gene',
                               right='ncbi_taxonomy_operon')
        return operon

    @staticmethod
    def transform_gene(operon: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        mask = operon['protrend_id'].notna()

        genes_in_db = operon[mask].copy()
        genes = operon[~mask].copy()

        genes = create_input_value(genes, 'locus_tag')
        return genes_in_db, genes

    def transform(self):
        conserved = read_from_stack(stack=self.transform_stack, key='conserved', columns=self.conserved_columns,
                                    reader=read_txt)
        known = read_from_stack(stack=self.transform_stack, key='known', columns=self.known_columns,
                                reader=read_txt)

        operon = self.transform_operon(conserved, known)
        gene = self.fetch_gene()
        operon = self.merge_operon_gene(operon, gene)

        genes_in_db, genes = self.transform_gene(operon)
        annotated_genes = self.annotate_genes(genes)

        df = pd.merge(annotated_genes, genes, on='input_value', suffixes=('_annotation', '_operondb'))

        # merge loci
        df = merge_columns(df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_operondb')

        df = df.drop(columns=['input_value', 'protrend_id'])

        df = pd.concat([df, genes_in_db])
        df = df.reset_index(drop=True)

        self.stack_transformed_nodes(df)
        return df


class GeneToOrganismConnector(OperonDBConnector,
                              source='operondb',
                              version='0.0.0',
                              from_node=Gene,
                              to_node=Organism,
                              register=True):
    default_connect_stack = {'gene': 'integrated_gene.json'}

    def connect(self):
        df = self.create_connection(source='gene', target='gene',
                                    target_column='organism')
        self.stack_connections(df)
