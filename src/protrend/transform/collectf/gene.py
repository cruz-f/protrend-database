import pandas as pd

from protrend.io import read_from_stack, read_json_lines, read_json_frame
from protrend.model import Gene
from protrend.transform.collectf.base import CollectfTransformer
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.utils import SetList
from protrend.utils.processors import (to_int_str, apply_processors, rstrip, lstrip, to_list_nan)


class GeneTransformer(CollectfTransformer,
                      source='collectf',
                      version='0.0.1',
                      node=Gene,
                      order=80,
                      register=True):
    default_transform_stack = {'gene': 'Gene.json', 'regulator': 'integrated_regulator.json'}
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop',
                       'regulon', 'operon', 'tfbs',
                       'regulator_uniprot_accession', 'ncbi_taxonomy', 'organism_protrend_id',
                       'locus_tag_old'])
    read_columns = SetList(['locus_tag', 'regulon', 'operon', 'tfbs'])

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = self.select_columns(regulator, 'uniprot_accession', 'ncbi_taxonomy', 'organism_protrend_id')
        regulator = regulator.rename(columns={'uniprot_accession': 'regulator_uniprot_accession'})
        regulator = apply_processors(regulator, ncbi_taxonomy=to_int_str)
        return regulator

    def transform_gene(self, gene: pd.DataFrame, regulator: pd.DataFrame) -> pd.DataFrame:
        gene = apply_processors(gene, locus_tag=[rstrip, lstrip], regulon=to_list_nan)
        gene = gene.dropna(subset=['locus_tag'])
        gene = self.drop_empty_string(gene, 'locus_tag')
        gene = self.drop_duplicates(df=gene, subset=['locus_tag'])

        gene = gene.explode('regulon')

        gene = pd.merge(gene, regulator, left_on='regulon', right_on='regulator_uniprot_accession')
        gene = self.drop_duplicates(df=gene, subset=['locus_tag', 'organism_protrend_id'], perfect_match=True)

        gene = gene.assign(locus_tag_old=gene['locus_tag'].copy())
        gene = self.create_input_value(df=gene, col='locus_tag')
        return gene

    def transform(self):
        # noinspection DuplicatedCode
        gene = read_from_stack(stack=self.transform_stack, key='gene',
                               columns=self.read_columns, reader=read_json_lines)

        regulator = read_from_stack(stack=self.transform_stack, key='regulator',
                                    columns=RegulatorTransformer.columns, reader=read_json_frame)

        regulator = self.transform_regulator(regulator)
        genes = self.transform_gene(gene=gene, regulator=regulator)
        annotated_genes = self.annotate_genes(genes)

        df = pd.merge(annotated_genes, genes, on='input_value', suffixes=('_annotation', '_collectf'))

        # merge loci
        df = self.merge_loci(df=df, left_suffix='_annotation', right_suffix='_collectf')

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
