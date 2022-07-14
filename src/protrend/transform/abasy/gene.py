import pandas as pd

from protrend.model import Gene
from protrend.transform.abasy.base import AbasyTransformer, read_abasy_genes
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.transformations import (create_input_value, drop_duplicates, drop_empty_string, merge_loci,
                                                merge_columns)
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_int_str


class GeneTransformer(GeneMixIn, AbasyTransformer,
                      source='abasy',
                      version='0.0.0',
                      node=Gene,
                      order=100,
                      register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'protein_sequence', 'strand', 'start', 'stop',
                       'Gene_name', 'Locus_tag', 'NCBI_gene_ID', 'Uniprot_ID', 'Synonyms',
                       'Product_function', 'NDA_component', 'taxonomy', 'source', 'ncbi_taxonomy', 'gene_taxonomy'])

    @staticmethod
    def transform_gene(gene: pd.DataFrame) -> pd.DataFrame:
        gene = gene.dropna(subset=['Gene_name', 'taxonomy'])
        gene = drop_empty_string(gene, 'Gene_name', 'taxonomy')
        gene = drop_duplicates(df=gene, subset=['Gene_name', 'taxonomy'], perfect_match=True)

        gene = apply_processors(gene,
                                Gene_name=[rstrip, lstrip],
                                Locus_tag=[rstrip, lstrip],
                                NCBI_gene_ID=[to_int_str, rstrip, lstrip],
                                Uniprot_ID=[rstrip, lstrip])

        gene_taxonomy = gene['Gene_name'] + gene['taxonomy']

        gene = gene.assign(gene_taxonomy=gene_taxonomy,
                           ncbi_taxonomy=gene['taxonomy'].copy(),
                           locus_tag=gene['Locus_tag'].copy(),
                           name=gene['Gene_name'].copy(),
                           ncbi_gene=gene['NCBI_gene_ID'].copy(),
                           uniprot_accession=gene['Uniprot_ID'].copy())

        gene = create_input_value(df=gene, col='gene_taxonomy')
        return gene

    def transform(self):
        gene = read_abasy_genes(self.source, self.version)

        genes = self.transform_gene(gene)
        annotated_genes = self.annotate_genes(genes)

        df = pd.merge(annotated_genes, genes, on='input_value', suffixes=('_annotation', '_abasy'))

        # merge loci
        df = merge_loci(df=df, left_suffix='_annotation', right_suffix='_abasy')

        # merge name
        df = merge_columns(df=df, column='name', left='name_annotation', right='name_abasy')

        # merge ncbi_gene
        df = merge_columns(df=df, column='ncbi_gene', left='ncbi_gene_annotation', right='ncbi_gene_abasy')
        pseudo_cap_mask = df['ncbi_gene'].str.startswith('PseudoCap')
        pseudo_cap_mask = pseudo_cap_mask.fillna(False)
        df.loc[pseudo_cap_mask, 'ncbi_gene'] = None

        # merge uniprot_accession
        df = merge_columns(df=df, column='uniprot_accession', left='uniprot_accession_annotation',
                           right='uniprot_accession_abasy')

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
