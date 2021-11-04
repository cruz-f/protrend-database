from typing import List, Union

import pandas as pd

from protrend.model.model import Gene
from protrend.transform import GeneDTO
from protrend.transform.abasy.base import AbasyTransformer
from protrend.transform.annotation import annotate_genes
from protrend.transform.processors import apply_processors, rstrip, lstrip, to_int_str
from protrend.utils import SetList


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
                       'Product_function', 'NDA_component', 'taxonomy', 'gene_name_taxonomy'])

    def _transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = self.drop_duplicates(df=gene, subset=['Gene_name', 'taxonomy'],
                                    perfect_match=True)
        gene = gene.dropna(subset=['Gene_name', 'taxonomy'])

        gene['gene_name_taxonomy'] = gene['Gene_name'] + gene['taxonomy']

        gene = apply_processors(gene,
                                gene_name_taxonomy=[rstrip, lstrip],
                                Gene_name=[rstrip, lstrip],
                                Locus_tag=[rstrip, lstrip],
                                NCBI_gene_ID=[to_int_str, rstrip, lstrip],
                                Uniprot_ID=[rstrip, lstrip])

        gene['locus_tag'] = gene['Locus_tag']
        gene['name'] = gene['Gene_name']
        gene['ncbi_gene'] = gene['NCBI_gene_ID']
        gene['uniprot_accession'] = gene['Uniprot_ID']

        gene = self.create_input_value(df=gene, col='gene_name_taxonomy')
        return gene

    @staticmethod
    def _annotate_genes(input_values: List[Union[None, str]],
                        loci: List[Union[None, str]],
                        names: List[str],
                        ncbi_genes: List[str],
                        uniprot_accessions: List[str],
                        taxa: List[str]):

        dtos = [GeneDTO(input_value=input_value) for input_value in input_values]
        annotate_genes(dtos=dtos, loci=loci, names=names, taxa=taxa,
                       uniprot_proteins=uniprot_accessions, ncbi_genes=ncbi_genes)

        for dto, name in zip(dtos, names):
            dto.synonyms.append(name)

        # locus_tag: List[str]
        # name: List[str]
        # synonyms: List[str]
        # function: List[str]
        # description: List[str]
        # ncbi_gene: List[str]
        # ncbi_protein: List[str]
        # genbank_accession: List[str]
        # refseq_accession: List[str]
        # uniprot_accession: List[str]
        # sequence: List[str]
        # strand: List[str]
        # start: List[int]
        # stop: List[int]

        genes = pd.DataFrame([dto.to_dict() for dto in dtos])
        strand_mask = (genes['strand'] != 'reverse') & (genes['strand'] != 'forward')
        genes.loc[strand_mask, 'strand'] = None
        return genes

    def transform(self):
        genes = self._build_genes()
        gene = self._transform_gene(genes)

        input_values = gene['input_value'].tolist()
        loci = gene['locus_tag'].tolist()
        names = gene['name'].tolist()
        ncbi_genes = gene['ncbi_gene'].tolist()
        uniprot_accessions = gene['uniprot_accession'].tolist()
        taxa = gene['taxonomy'].tolist()

        annotated_genes = self._annotate_genes(input_values=input_values, loci=loci, names=names, taxa=taxa,
                                               ncbi_genes=ncbi_genes, uniprot_accessions=uniprot_accessions)

        df = pd.merge(annotated_genes, gene, on='input_value', suffixes=('_annotation', '_abasy'))

        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_abasy')
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_abasy')

        df = self.merge_columns(df=df, column='ncbi_gene', left='ncbi_gene_annotation', right='ncbi_gene_abasy')
        pseudo_cap_mask = df['ncbi_gene'].str.startswith('PseudoCap')
        pseudo_cap_mask = pseudo_cap_mask.fillna(False)
        df.loc[pseudo_cap_mask, 'ncbi_gene'] = None

        df = self.merge_columns(df=df, column='uniprot_accession', left='uniprot_accession_annotation',
                                right='uniprot_accession_abasy')

        df = df.drop(columns=['input_value'])

        self._stack_transformed_nodes(df)

        return df
