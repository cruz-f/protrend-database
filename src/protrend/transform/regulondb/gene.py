import pandas as pd
from Bio import SeqIO

from protrend.io import read_from_stack
from protrend.model.model import Gene
from protrend.transform.regulondb.base import RegulondbTransformer, regulondb_reader
from protrend.transform import transform_sequence
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip


class GeneTransformer(RegulondbTransformer,
                      source='regulondb',
                      version='0.0.0',
                      node=Gene,
                      order=100,
                      register=True):
    default_transform_stack = {'gene': 'gene.txt', 'sequence': 'sequence.gb'}
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop',
                       'name_lower',
                       'gene_id', 'gene_name', 'gene_posleft', 'gene_posright', 'gene_strand',
                       'gene_sequence', 'gc_content', 'cri_score', 'gene_note',
                       'gene_internal_comment', 'key_id_org', 'gene_type'])
    read_columns = SetList(['gene_id', 'gene_name', 'gene_posleft', 'gene_posright', 'gene_strand',
                            'gene_sequence', 'gc_content', 'cri_score', 'gene_note',
                            'gene_internal_comment', 'key_id_org', 'gene_type'])

    def transform_gene(self, gene: pd.DataFrame, sequence: pd.DataFrame) -> pd.DataFrame:
        gene = gene.assign(name=gene['gene_name'].copy())

        # noinspection DuplicatedCode
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
        reader = regulondb_reader(skiprows=39, names=self.read_columns)
        # noinspection DuplicatedCode
        gene = read_from_stack(stack=self.transform_stack, key='gene',
                               columns=self.read_columns, reader=reader)

        gb_file = self.transform_stack['sequence']
        sequence = SeqIO.read(gb_file, "genbank")
        sequence = transform_sequence(sequence)

        genes = self.transform_gene(gene=gene, sequence=sequence)
        annotated_genes = self.annotate_genes(genes)

        df = pd.merge(annotated_genes, genes, on='input_value', suffixes=('_annotation', '_regulondb'))

        df = self.merge_columns(df=df, column='locus_tag',
                                left='locus_tag_annotation', right='locus_tag_regulondb')
        df = self.merge_columns(df=df, column='name',
                                left='name_annotation', right='name_regulondb')
        df = self.merge_columns(df=df, column='genbank_accession',
                                left='genbank_accession_annotation', right='genbank_accession_regulondb')
        df = self.merge_columns(df=df, column='uniprot_accession',
                                left='uniprot_accession_annotation', right='uniprot_accession_regulondb')

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
