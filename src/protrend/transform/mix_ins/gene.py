import pandas as pd

from protrend.annotation import (GeneDTO, annotate_genes)
from protrend.log import ProtrendLogger
from ._utils import get_values
from protrend.utils.constants import UNKNOWN, REVERSE, FORWARD
from protrend.transform.transformations import drop_empty_string, merge_columns
from protrend.io import read_json_frame
from protrend.utils import Settings


class GeneMixIn:

    @staticmethod
    def annotate_genes(df: pd.DataFrame) -> pd.DataFrame:
        input_values = get_values(df, 'input_value')

        genes = [GeneDTO(input_value=input_value) for input_value in input_values]

        ProtrendLogger.log.info(f'Annotating {len(genes)} genes')

        loci = get_values(df, 'locus_tag')
        names = get_values(df, 'name')
        taxa = get_values(df, 'ncbi_taxonomy')
        uniprot_proteins = get_values(df, 'uniprot_accession')
        ncbi_proteins = get_values(df, 'ncbi_protein')
        ncbi_genbanks = get_values(df, 'genbank_accession')
        ncbi_refseqs = get_values(df, 'refseq_accession')
        ncbi_genes = get_values(df, 'ncbi_gene')

        iterator = zip(
            ('locus_tag', 'name', 'ncbi_taxonomy', 'uniprot_accession', 'ncbi_protein', 'genbank_accession',
             'refseq_accession', 'ncbi_gene'),
            (loci, names, taxa, uniprot_proteins, ncbi_proteins, ncbi_genbanks, ncbi_refseqs, ncbi_genes)
        )

        params = [param for param, value in iterator if value is not None]
        params = ','.join(params)

        ProtrendLogger.log.info(f'Annotating with the following params: {params}')

        annotate_genes(dtos=genes,
                       loci=loci,
                       names=names,
                       taxa=taxa,
                       uniprot_proteins=uniprot_proteins,
                       ncbi_proteins=ncbi_proteins,
                       ncbi_genbanks=ncbi_genbanks,
                       ncbi_refseqs=ncbi_refseqs,
                       ncbi_genes=ncbi_genes)

        genes_dict = [dto.to_dict() for dto in genes]
        genes_df = pd.DataFrame(genes_dict)

        if genes_df.empty:
            return genes_df

        genes_df = genes_df.dropna(subset=['locus_tag'])
        genes_df = drop_empty_string(genes_df, 'locus_tag')

        strand_mask = (genes_df['strand'] != REVERSE) & (genes_df['strand'] != FORWARD)
        genes_df.loc[strand_mask, 'strand'] = UNKNOWN

        return genes_df

    def annotate_with_genomes_database(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotates the given dataframe with the genomes database.
        """
        taxa = get_values(df, 'ncbi_taxonomy')

        annotated_genes = self.annotate_genes(df)

        ProtrendLogger.log.info('Annotating with the genomes database')

        genome_path = Settings.genomes_database.joinpath(f'{taxa}.json')

        if not genome_path.exists():
            return annotated_genes

        annotated_genes['synonyms'] = annotated_genes['synonyms'].fillna("").apply(list)
        largest_annotated_synonyms = annotated_genes['synonyms'].str.len().max()
        annotated_synonyms_cols = [f'synonym_{i+1}' for i in range(largest_annotated_synonyms)]
        genes_synonyms = pd.DataFrame(annotated_genes['synonyms'].to_list(), columns=annotated_synonyms_cols)
        annotated_genes = annotated_genes.drop(columns=['synonyms'])
        annotated_genes = annotated_genes.join(genes_synonyms)

        genome = read_json_frame(genome_path)
        genome['synonyms'] = genome['synonyms'].fillna("").apply(list)
        largest_genome_synonyms = genome['synonyms'].str.len().max()
        genome_synonyms_cols = [f'synonym_{i+1}' for i in range(largest_genome_synonyms)]
        genome_synonyms = pd.DataFrame(genome['synonyms'].to_list(), columns=genome_synonyms_cols)
        genome = genome.drop(columns=['synonyms'])
        genome = genome.join(genome_synonyms)

        left_on_cols = ['locus_tag', 'genbank_accession', 'uniprot_accession'] + annotated_synonyms_cols
        right_on_cols = ['locus_tag', 'genbank_accession', 'uniprot_accession'] + genome_synonyms_cols
        df = pd.merge(annotated_genes,
                      genome,
                      how='left',
                      left_on=left_on_cols,
                      right_on=right_on_cols,
                      suffixes=('_annotated', '_genomes'))

        df = df.drop(columns=['name_annotated',
                              'promoter_sequence',
                              'promoter_start',
                              'promoter_end',
                              'promoter_strand'])
        df = df.rename(columns={'name_genomes': 'name'})

        df = merge_columns(df, column='protein_sequence', left='protein_sequence', right='sequence')
        df = merge_columns(df, column='strand', left='gene_strand', right='strand')
        df = merge_columns(df, column='start', left='gene_start', right='start')
        df = merge_columns(df, column='stop', left='gene_stop', right='stop')
        return df
