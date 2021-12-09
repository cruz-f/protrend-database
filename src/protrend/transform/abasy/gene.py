import pandas as pd

from protrend.io import read_csv
from protrend.model import Gene
from protrend.transform.abasy.base import AbasyTransformer
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_int_str


class GeneTransformer(AbasyTransformer,
                      source='abasy',
                      version='0.0.0',
                      node=Gene,
                      order=100,
                      register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop',
                       'Gene_name', 'Locus_tag', 'NCBI_gene_ID', 'Uniprot_ID', 'Synonyms',
                       'Product_function', 'NDA_component', 'taxonomy', 'gene_taxonomy'])

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = self.drop_duplicates(df=gene, subset=['Gene_name', 'taxonomy'], perfect_match=True)
        gene = gene.dropna(subset=['Gene_name', 'taxonomy'])
        gene = self.drop_empty_string(gene, 'Gene_name')

        gene = apply_processors(gene,
                                Gene_name=[rstrip, lstrip],
                                Locus_tag=[rstrip, lstrip],
                                NCBI_gene_ID=[to_int_str, rstrip, lstrip],
                                Uniprot_ID=[rstrip, lstrip])

        gene_taxonomy = gene['Gene_name'] + gene['taxonomy']

        gene = gene.assign(gene_taxonomy=gene_taxonomy,
                           locus_tag=gene['Locus_tag'].copy(),
                           name=gene['Gene_name'].copy(),
                           ncbi_gene=gene['NCBI_gene_ID'].copy(),
                           uniprot_accession=gene['Uniprot_ID'].copy())

        gene = self.create_input_value(df=gene, col='gene_taxonomy')
        return gene

    def transform(self):
        gene_stack = self.contact_stacks(stack=self.gene_stack,
                                         taxa=self.taxa_to_organism_code,
                                         default_columns=self.default_gene_columns,
                                         reader=read_csv)

        genes = self.transform_gene(gene_stack)
        annotated_genes = self.annotate_genes(genes)

        df = pd.merge(annotated_genes, genes, on='input_value', suffixes=('_annotation', '_abasy'))

        # merge loci
        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_abasy')
        df = df.dropna(subset=['locus_tag'])
        df = self.drop_empty_string(df, 'locus_tag')
        df = self.drop_duplicates(df=df, subset=['locus_tag'], perfect_match=True)

        # merge name
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_abasy')

        # merge ncbi_gene
        df = self.merge_columns(df=df, column='ncbi_gene', left='ncbi_gene_annotation', right='ncbi_gene_abasy')
        pseudo_cap_mask = df['ncbi_gene'].str.startswith('PseudoCap')
        pseudo_cap_mask = pseudo_cap_mask.fillna(False)
        df.loc[pseudo_cap_mask, 'ncbi_gene'] = None

        # merge uniprot_accession
        df = self.merge_columns(df=df, column='uniprot_accession', left='uniprot_accession_annotation',
                                right='uniprot_accession_abasy')

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)

        return df
