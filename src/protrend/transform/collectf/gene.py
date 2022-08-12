import pandas as pd

from protrend.io.utils import read_regulator
from protrend.model import Gene
from protrend.report import ProtrendReporter
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
                       'regulon_id', 'url', 'operon', 'gene', 'tfbs', 'experimental_evidence',
                       'organism_protrend_id', 'organism_name', 'ncbi_taxonomy'])

    @staticmethod
    def transform_regulator(regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = select_columns(regulator, 'regulon_id', 'url', 'operon', 'gene', 'tfbs', 'experimental_evidence',
                                   'organism_protrend_id', 'organism_name', 'ncbi_taxonomy')
        regulator = apply_processors(regulator, ncbi_taxonomy=to_int_str)
        return regulator

    @staticmethod
    def transform_gene(regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = apply_processors(regulator, gene=to_list_nan)

        gene = regulator.explode('gene')

        gene = apply_processors(gene, gene=[rstrip, lstrip])
        gene = gene.dropna(subset=['gene'])
        gene = drop_empty_string(gene, 'gene')
        gene = drop_duplicates(df=gene, subset=['gene'])

        gene = gene.assign(locus_tag=gene['gene'].copy())

        gene = create_input_value(df=gene, col='locus_tag')
        return gene

    def transform(self):
        regulator = read_regulator(source=self.source, version=self.version, columns=RegulatorTransformer.columns)

        regulator = self.transform_regulator(regulator)
        genes = self.transform_gene(regulator=regulator)

        ProtrendReporter.report_objects(source=self.source, version=self.version,
                                        system='extract', label=self.node.node_name(),
                                        objects=genes.shape[0], properties=genes.shape[1])

        annotated_genes = self.annotate_genes(genes)

        df = pd.merge(annotated_genes, genes, on='input_value', suffixes=('_annotation', '_collectf'))

        # merge loci
        df = merge_loci(df=df, left_suffix='_annotation', right_suffix='_collectf')

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
