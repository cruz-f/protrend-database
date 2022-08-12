import pandas as pd

from protrend.io import read_json_lines, read
from protrend.io.utils import read_regulator
from protrend.model import Gene
from protrend.report import ProtrendReporter
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.regprecise.base import RegPreciseTransformer
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.transformations import (select_columns, drop_empty_string, group_by, create_input_value,
                                                merge_columns)
from protrend.utils import SetList
from protrend.utils.processors import (rstrip, lstrip, apply_processors, take_last,
                                       flatten_set_list_nan, to_int_str, to_list_nan, to_set_list)


class GeneTransformer(GeneMixIn, RegPreciseTransformer,
                      source='regprecise',
                      version='0.0.0',
                      node=Gene,
                      order=80,
                      register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'protein_sequence', 'strand', 'start', 'stop',
                       'url', 'regulon', 'operon', 'tfbs',
                       'ncbi_taxonomy', 'regprecise_locus_tag'])

    @staticmethod
    def transform_regulator(regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = select_columns(regulator, 'regulon_id', 'ncbi_taxonomy')
        regulator = regulator.rename(columns={'regulon_id': 'regulon'})
        regulator = apply_processors(regulator, regulon=to_int_str, ncbi_taxonomy=to_int_str)
        return regulator

    @staticmethod
    def transform_gene(gene: pd.DataFrame, regulator: pd.DataFrame) -> pd.DataFrame:
        gene = gene.assign(regprecise_locus_tag=gene['locus_tag'].copy())

        gene = gene.dropna(subset=['locus_tag'])
        gene = drop_empty_string(gene, 'locus_tag')

        gene = apply_processors(gene,
                                locus_tag=[rstrip, lstrip], name=[rstrip, lstrip], function=[rstrip, lstrip],
                                regulon=to_list_nan)

        aggregation = {'name': take_last, 'function': take_last, 'url': to_set_list,
                       'regprecise_locus_tag': to_set_list}
        gene = group_by(df=gene, column='locus_tag', aggregation=aggregation, default=flatten_set_list_nan)

        gene = gene.explode('regulon')
        gene = apply_processors(gene, regulon=to_int_str)

        # + 'regulon', 'ncbi_taxonomy'
        gene = pd.merge(gene, regulator, on='regulon')

        aggregation = {'name': take_last, 'function': take_last, 'ncbi_taxonomy': take_last, 'regulon': to_set_list}
        gene = group_by(df=gene, column='locus_tag', aggregation=aggregation, default=flatten_set_list_nan)

        gene = create_input_value(df=gene, col='locus_tag')
        return gene

    def transform(self):
        # noinspection DuplicatedCode
        gene = read(source=self.source, version=self.version,
                    file='Gene.json', reader=read_json_lines,
                    default=pd.DataFrame(columns=['locus_tag', 'name', 'function', 'url',
                                                  'regulon', 'operon', 'tfbs']))

        regulator = read_regulator(source=self.source, version=self.version, columns=RegulatorTransformer.columns)

        regulator = self.transform_regulator(regulator)

        ProtrendReporter.report_objects(source=self.source, version=self.version,
                                        system='extract', label=self.node.node_name(),
                                        objects=regulator.shape[0], properties=regulator.shape[1])

        genes = self.transform_gene(gene=gene, regulator=regulator)

        annotated_genes = self.annotate_genes(genes)

        df = pd.merge(annotated_genes, genes, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_regprecise')
        df = merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')
        df = merge_columns(df=df, column='function', left='function_annotation', right='function_regprecise')

        df = apply_processors(df, regulon=to_int_str, ncbi_taxonomy=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
