import pandas as pd

from protrend.io import read_json_lines, read_json_frame, read_from_stack
from protrend.model import Gene
from protrend.transform.regprecise.base import RegPreciseTransformer
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.utils import SetList
from protrend.utils.processors import (rstrip, lstrip, apply_processors, take_last,
                                       flatten_set_list, to_int_str, to_list_nan, to_set_list)


class GeneTransformer(RegPreciseTransformer,
                      source='regprecise',
                      version='0.0.0',
                      node=Gene,
                      order=80,
                      register=True):
    default_transform_stack = {'gene': 'Gene.json', 'regulator': 'integrated_regulator.json'}
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop',
                       'url', 'regulon', 'operon', 'tfbs',
                       'ncbi_taxonomy', 'regprecise_locus_tag'])
    read_columns = SetList(['locus_tag', 'name', 'function', 'url', 'regulon', 'operon', 'tfbs'])

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = self.select_columns(regulator, 'regulon_id', 'ncbi_taxonomy')
        regulator = regulator.rename(columns={'regulon_id': 'regulon'})
        regulator = apply_processors(regulator, regulon=to_int_str, ncbi_taxonomy=to_int_str)
        return regulator

    def transform_gene(self, gene: pd.DataFrame, regulator: pd.DataFrame) -> pd.DataFrame:
        gene = gene.assign(regprecise_locus_tag=gene['locus_tag'].copy())

        gene = gene.dropna(subset=['locus_tag'])
        gene = self.drop_empty_string(gene, 'locus_tag')

        gene = apply_processors(gene,
                                locus_tag=[rstrip, lstrip], name=[rstrip, lstrip], function=[rstrip, lstrip],
                                regulon=to_list_nan)

        aggregation = {'name': take_last, 'function': take_last, 'url': to_set_list,
                       'regprecise_locus_tag': to_set_list}
        gene = self.group_by(df=gene, column='locus_tag', aggregation=aggregation, default=flatten_set_list)

        gene = gene.explode('regulon')
        gene = apply_processors(gene, regulon=to_int_str)

        # + 'regulon', 'ncbi_taxonomy'
        gene = pd.merge(gene, regulator, on='regulon')

        aggregation = {'name': take_last, 'function': take_last, 'ncbi_taxonomy': take_last, 'regulon': to_set_list}
        gene = self.group_by(df=gene, column='locus_tag', aggregation=aggregation, default=flatten_set_list)

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

        df = pd.merge(annotated_genes, genes, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_regprecise')
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')
        df = self.merge_columns(df=df, column='function', left='function_annotation', right='function_regprecise')

        df = apply_processors(df, regulon=to_int_str, ncbi_taxonomy=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
