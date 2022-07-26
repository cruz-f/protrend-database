from typing import Optional, Tuple

import pandas as pd

from protrend.annotation import (GeneDTO, annotate_genes)
from protrend.io import read_json_frame
from protrend.log import ProtrendLogger
from protrend.transform.transformations import drop_empty_string, merge_columns, group_by
from protrend.utils import Settings
from protrend.utils.constants import UNKNOWN, REVERSE, FORWARD
from ._utils import get_values
from ...utils.processors import take_first, to_set_list, take_last, flatten_set_list_nan


def _explode_synonyms(df: pd.DataFrame) -> pd.DataFrame:
    df['synonyms'] = df['synonyms'].fillna("").apply(list)

    df = df.explode(column='synonyms')
    df = df.assign(synonyms=df['synonyms'].str.lower())
    df = df.drop_duplicates(subset=['synonyms'])
    return df


class GeneMixIn:

    @staticmethod
    def _read_genome(taxa: str) -> Optional[pd.DataFrame]:
        """
        Reads the genomes database for the given taxa.
        """
        genome_path = Settings.genomes_database.joinpath(f'{taxa}.json')

        if not genome_path.exists():
            ProtrendLogger.log.info(f'Missing genome database for taxonomy {taxa}')
            return pd.DataFrame(columns=['locus_tag', 'name', 'synonyms', 'uniprot_accession',
                                         'genbank_accession', 'gene_sequence', 'start',
                                         'end', 'strand'])

        genome = read_json_frame(genome_path)
        genome = genome.drop(columns=['promoter_sequence',
                                      'promoter_start',
                                      'promoter_end',
                                      'promoter_strand'])
        genome = genome.rename(columns={'gene_strand': 'strand',
                                        'gene_start': 'start',
                                        'gene_end': 'stop'})
        genome = genome.assign(strand=genome['strand'].map({1: FORWARD, -1: REVERSE}))
        genome = genome.dropna().drop_duplicates(subset=['locus_tag']).reset_index(drop=True)

        ProtrendLogger.log.info(f'Loaded genome database for taxonomy {taxa}')
        return genome

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

    def _annotate_with_genomes_database(self, annotated_genes: pd.DataFrame, taxa: str) -> pd.DataFrame:
        """
        Annotates the given dataframe with the genomes database.
        """
        genome = self._read_genome(taxa)
        if genome.empty:
            return annotated_genes

        genome = _explode_synonyms(genome)

        annotated_genes = annotated_genes.reset_index(drop=True)
        annotated_genes = _explode_synonyms(annotated_genes)

        annotated_genes_cols = list(annotated_genes.columns)

        for col in ['synonyms', 'uniprot_accession', 'locus_tag']:

            annotated_genes = annotated_genes.merge(genome, how='left', on=col,
                                                    suffixes=('_annotated_x', '_genome_y'))

            for merged_col in annotated_genes_cols:
                left_col = f'{merged_col}_annotated_x'
                right_col = f'{merged_col}_genome_y'

                if left_col in annotated_genes.columns and right_col in annotated_genes.columns:
                    # we want to keep the genome annotation rather than the annotation obtained
                    # by searching the NCBI and UniProt
                    annotated_genes = merge_columns(annotated_genes, column=merged_col,
                                                    right=left_col, left=right_col)

        annotated_genes = annotated_genes.reset_index(drop=True)

        aggregation = {
            'synonyms': flatten_set_list_nan,
        }

        annotated_genes = group_by(annotated_genes,
                                   column='input_value',
                                   aggregation=aggregation,
                                   default=take_last)
        annotated_genes = annotated_genes.reset_index(drop=True)
        return annotated_genes

    def annotate_genes(self, df: pd.DataFrame, use_genomes_database: bool = True) -> pd.DataFrame:
        """
        Annotates the given dataframe using NCBI and Uniprot databases.
        """
        annotated_genes = self._annotate_genes(df)

        if use_genomes_database:
            unique_taxa = annotated_genes['ncbi_taxonomy'].unique()
            dfs = []
            for taxa in unique_taxa:
                annotated_genes_by_taxa = annotated_genes[annotated_genes['ncbi_taxonomy'] == taxa]
                annotated_genes_by_taxa = annotated_genes_by_taxa.reset_index(drop=True)

                df = self._annotate_with_genomes_database(annotated_genes_by_taxa, taxa)
                dfs.append(df)

            annotated_genes = pd.concat(dfs, ignore_index=True)
            annotated_genes = annotated_genes.reset_index(drop=True)

        annotated_genes = annotated_genes.drop(columns=['ncbi_taxonomy'])
        return annotated_genes
