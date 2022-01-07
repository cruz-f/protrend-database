import pandas as pd
from Bio import SeqIO

from protrend.io import read_json_lines, read_from_stack
from protrend.model import Gene
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.transform.transformer import transform_sequence
from protrend.utils import SetList
from protrend.utils.processors import (rstrip, lstrip, apply_processors)


class GeneTransformer(DBTBSTransformer,
                      source='dbtbs',
                      version='0.0.4',
                      node=Gene,
                      order=100,
                      register=True):
    default_transform_stack = {'gene': 'Gene.json', 'sequence': 'sequence.gb'}
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop',
                       'url', 'regulation', 'pubmed', 'tf', 'tfbs',
                       'name_lower', 'name_dbtbs'])
    read_columns = SetList(['name', 'url', 'regulation', 'pubmed', 'tf', 'tfbs'])

    def transform_gene(self, gene: pd.DataFrame, sequence: pd.DataFrame) -> pd.DataFrame:
        gene = gene.explode(column='name')
        gene = gene.explode(column='url')
        gene = gene.explode(column='regulation')
        gene = gene.explode(column='tf')
        gene = gene.explode(column='tfbs')

        gene = gene.assign(name_dbtbs=gene['name'].copy())

        gene = apply_processors(gene, name=[rstrip, lstrip])

        # filter nan and duplicates
        gene = gene.dropna(subset=['name'])
        gene = self.drop_empty_string(gene, 'name')
        gene = self.drop_duplicates(df=gene, subset=['name'])

        gene = gene.assign(name_lower=gene['name'].str.lower())

        gene = pd.merge(gene, sequence, on='name_lower')

        gene = self.create_input_value(df=gene, col='locus_tag')
        return gene

    def transform(self):
        # noinspection DuplicatedCode
        gene = read_from_stack(stack=self.transform_stack, key='gene',
                               columns=self.read_columns, reader=read_json_lines)

        gb_file = self.transform_stack['sequence']
        sequence = SeqIO.read(gb_file, "genbank")
        sequence = transform_sequence(sequence)

        genes = self.transform_gene(gene=gene, sequence=sequence)
        annotated_genes = self.annotate_genes(genes)

        df = self.merge_annotations(annotated_genes, genes)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
