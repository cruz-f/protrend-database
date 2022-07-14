import pandas as pd

from protrend.io import read, read_genbank
from protrend.model import Gene
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.regulondb.base import RegulonDBTransformer, regulondb_reader
from protrend.transform.transformations import drop_empty_string, drop_duplicates, create_input_value, merge_columns
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip


class GeneTransformer(GeneMixIn, RegulonDBTransformer,
                      source='regulondb',
                      version='0.0.0',
                      node=Gene,
                      order=100,
                      register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'protein_sequence', 'strand', 'start', 'stop',
                       'name_lower',
                       'gene_id', 'gene_name', 'gene_posleft', 'gene_posright', 'gene_strand',
                       'gene_sequence', 'gc_content', 'cri_score', 'gene_note',
                       'gene_internal_comment', 'key_id_org', 'gene_type'])

    @staticmethod
    def transform_gene(gene: pd.DataFrame, sequence: pd.DataFrame) -> pd.DataFrame:
        gene = gene.assign(name=gene['gene_name'].copy())

        # noinspection DuplicatedCode
        gene = apply_processors(gene, name=[rstrip, lstrip])

        # filter nan and duplicates
        gene = gene.dropna(subset=['name'])
        gene = drop_empty_string(gene, 'name')
        gene = drop_duplicates(df=gene, subset=['name'])

        gene = gene.assign(name_lower=gene['name'].str.lower())

        gene = pd.merge(gene, sequence, on='name_lower')

        gene = create_input_value(df=gene, col='locus_tag')
        return gene

    def transform(self):
        columns = ['gene_id', 'gene_name', 'gene_posleft', 'gene_posright', 'gene_strand',
                   'gene_sequence', 'gc_content', 'cri_score', 'gene_note',
                   'gene_internal_comment', 'key_id_org', 'gene_type']
        reader = regulondb_reader(skiprows=39, names=columns)
        gene = read(source=self.source, version=self.version, file='gene.txt', reader=reader,
                    default=pd.DataFrame(columns=columns))

        sequence = read(source=self.source, version=self.version, file='sequence.gb', reader=read_genbank,
                        default=pd.DataFrame(columns=['name_lower', 'locus_tag', 'genbank_accession',
                                                      'uniprot_accession']))

        genes = self.transform_gene(gene=gene, sequence=sequence)
        annotated_genes = self.annotate_genes(genes)

        df = pd.merge(annotated_genes, genes, on='input_value', suffixes=('_annotation', '_regulondb'))

        df = merge_columns(df=df, column='locus_tag',
                           left='locus_tag_annotation', right='locus_tag_regulondb')
        df = merge_columns(df=df, column='name',
                           left='name_annotation', right='name_regulondb')
        df = merge_columns(df=df, column='genbank_accession',
                           left='genbank_accession_annotation', right='genbank_accession_regulondb')
        df = merge_columns(df=df, column='uniprot_accession',
                           left='uniprot_accession_annotation', right='uniprot_accession_regulondb')

        df = df.dropna(subset=['locus_tag'])
        df = drop_empty_string(df, 'locus_tag')
        df = drop_duplicates(df, subset=['locus_tag'])

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
