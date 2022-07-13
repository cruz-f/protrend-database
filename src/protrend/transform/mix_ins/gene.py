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
    def _annotate_genes(df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotates the given dataframe with the genomes database.
        The output contains the ncbi taxonomy column!!!
        """
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
            genes_df = genes_df.assign(ncbi_taxonomy=None)
            return genes_df

        genes_df = genes_df.assign(ncbi_taxonomy=taxa)
        genes_df = genes_df.dropna(subset=['locus_tag'])
        genes_df = drop_empty_string(genes_df, 'locus_tag')

        strand_mask = (genes_df['strand'] != REVERSE) & (genes_df['strand'] != FORWARD)
        genes_df.loc[strand_mask, 'strand'] = UNKNOWN

        return genes_df

    def annotate_genes(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotates the given dataframe using NCBI and Uniprot databases.
        """
        genes_df = self._annotate_genes(df)
        genes_df = genes_df.drop(columns=['ncbi_taxonomy'])
        return genes_df

    @staticmethod
    def _annotate_with_genomes_database(annotated_genes: pd.DataFrame, taxa: str) -> pd.DataFrame:
        """
        Annotates the given dataframe with the genomes database.
        """
        genome_path = Settings.genomes_database.joinpath(f'{taxa}.json')

        if not genome_path.exists():
            return annotated_genes

        # expanding the df so that each synonym is a separate column
        annotated_genes['synonyms'] = annotated_genes['synonyms'].fillna("").apply(list)
        largest_synonyms_list = annotated_genes['synonyms'].str.len().max()
        synonyms_cols = [f'synonym_{i + 1}' for i in range(largest_synonyms_list)]
        synonyms = pd.DataFrame(annotated_genes['synonyms'].to_list(), columns=synonyms_cols)

        annotated_genes = annotated_genes.drop(columns=['synonyms']).join(synonyms)

        genome = read_json_frame(genome_path)
        genome = genome.drop(columns=['promoter_sequence',
                                      'promoter_start',
                                      'promoter_end',
                                      'promoter_strand'])
        genome = genome.rename(columns={'gene_strand': 'strand',
                                        'gene_start': 'start',
                                        'gene_end': 'stop'})
        genome = genome.assign(strand=genome['strand'].map({1: FORWARD, -1: REVERSE}))

        # Repeating the locus tag to match the gene synonyms obtained during the annotation
        genome_synonyms = pd.concat([genome['locus_tag']] * largest_synonyms_list, axis=1, ignore_index=True)
        genome_synonyms.columns = [f'{synonyms_col}_locus_tag' for synonyms_col in synonyms_cols]
        genome = genome.join(genome_synonyms)

        df = annotated_genes.copy()
        for col in ['locus_tag', 'uniprot_accession']:
            df = df.set_index(col, drop=False)
            mapper = genome.dropna(subset=[col]).drop_duplicates(subset=[col]).set_index(col, drop=False)
            df.update(mapper)

        df = df.reset_index(drop=True)

        for col in synonyms_cols:
            df = df.set_index(col, drop=False)
            mapper = genome.dropna(subset=[f'{col}_locus_tag']).drop_duplicates(subset=[f'{col}_locus_tag']).\
                set_index([f'{col}_locus_tag'], drop=False)
            df.update(mapper)

        df = df.reset_index(drop=True)
        df = df.assign(synonyms=df[synonyms_cols].values.tolist())
        df['synonyms'] = df['synonyms'].apply(lambda x: set([y for y in x if y]))
        df = df.drop(columns=synonyms_cols)
        return df

    def annotate_with_genomes_database(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotates the given dataframe with the genomes database.
        """
        annotated_genes = self._annotate_genes(df)

        unique_taxa = annotated_genes['ncbi_taxonomy'].unique()
        dfs = []
        for taxa in unique_taxa:
            annotated_genes_by_taxa = annotated_genes[annotated_genes['ncbi_taxonomy'] == taxa]
            annotated_genes_by_taxa = annotated_genes_by_taxa.drop(columns=['ncbi_taxonomy']).reset_index(drop=True)
            df = self._annotate_with_genomes_database(annotated_genes_by_taxa, taxa)
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)
        df = df.reset_index(drop=True)
        return df
