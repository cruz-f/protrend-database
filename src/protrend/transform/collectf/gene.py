import pandas as pd

from protrend.io import read_json_lines, read
from protrend.io.utils import read_regulator
from protrend.model import Gene
from protrend.transform.collectf.base import CollecTFTransformer
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.transformations import (drop_empty_string, drop_duplicates, select_columns, create_input_value,
                                                merge_loci)
from protrend.utils import SetList
from protrend.utils.processors import (to_int_str, apply_processors, rstrip, lstrip, to_list_nan)


class GeneTransformer(GeneMixIn, CollecTFTransformer,
                      source='collectf',
                      version='0.0.1',
                      node=Gene,
                      order=80,
                      register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'protein_sequence', 'strand', 'start', 'stop',
                       'regulon', 'operon', 'tfbs',
                       'regulator_uniprot_accession', 'ncbi_taxonomy', 'organism_protrend_id',
                       'locus_tag_old'])

    @staticmethod
    def transform_regulator(regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = select_columns(regulator, 'uniprot_accession', 'ncbi_taxonomy', 'organism_protrend_id')
        regulator = regulator.rename(columns={'uniprot_accession': 'regulator_uniprot_accession'})
        regulator = apply_processors(regulator, ncbi_taxonomy=to_int_str)
        return regulator

    @staticmethod
    def transform_gene(gene: pd.DataFrame, regulator: pd.DataFrame) -> pd.DataFrame:
        gene = apply_processors(gene, locus_tag=[rstrip, lstrip], regulon=to_list_nan)
        gene = gene.dropna(subset=['locus_tag'])
        gene = drop_empty_string(gene, 'locus_tag')
        gene = drop_duplicates(df=gene, subset=['locus_tag'])

        gene = gene.explode('regulon')

        gene = pd.merge(gene, regulator, left_on='regulon', right_on='regulator_uniprot_accession')
        gene = drop_duplicates(df=gene, subset=['locus_tag', 'organism_protrend_id'], perfect_match=True)

        gene = gene.assign(locus_tag_old=gene['locus_tag'].copy())
        gene = create_input_value(df=gene, col='locus_tag')
        return gene

    def transform(self):
        gene = read(source=self.source, version=self.version,
                    file='Gene.json', reader=read_json_lines,
                    default=pd.DataFrame(columns=['locus_tag', 'regulon', 'operon', 'tfbs']))

        regulator = read_regulator(source=self.source, version=self.version, columns=RegulatorTransformer.columns)

        regulator = self.transform_regulator(regulator)
        genes = self.transform_gene(gene=gene, regulator=regulator)
        annotated_genes = self.annotate_genes(genes)

        df = pd.merge(annotated_genes, genes, on='input_value', suffixes=('_annotation', '_collectf'))

        # merge loci
        df = merge_loci(df=df, left_suffix='_annotation', right_suffix='_collectf')

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
