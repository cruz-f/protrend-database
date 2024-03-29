import pandas as pd

from protrend.io import read_json_lines, read, read_json_frame
from protrend.model import Gene
from protrend.report import ProtrendReporter
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.transformations import drop_empty_string, drop_duplicates, create_input_value
from protrend.utils import SetList, Settings
from protrend.utils.processors import rstrip, lstrip, apply_processors


class GeneTransformer(GeneMixIn, DBTBSTransformer,
                      source='dbtbs',
                      version='0.0.4',
                      node=Gene,
                      order=100,
                      register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'protein_sequence', 'strand', 'start', 'stop',
                       'url', 'regulation', 'pubmed', 'tf', 'tfbs',
                       'dbtbs_name'])

    @staticmethod
    def transform_gene(gene: pd.DataFrame, genome: pd.DataFrame) -> pd.DataFrame:
        gene = gene.explode(column='name')
        gene = gene.explode(column='url')
        gene = gene.explode(column='regulation')
        gene = gene.explode(column='tf')
        gene = gene.explode(column='tfbs')

        gene = gene.assign(dbtbs_name=gene['name'].copy())

        gene = apply_processors(gene, name=[rstrip, lstrip])

        # filter nan and duplicates
        gene = gene.dropna(subset=['name'])
        gene = drop_empty_string(gene, 'name')
        gene = drop_duplicates(df=gene, subset=['name'])

        gene = gene.assign(name_lower=gene['name'].str.lower())

        gene = pd.merge(gene, genome, on='name_lower')
        gene = gene.drop(columns=['name_lower'])

        # for locus tag annotation
        gene = gene.assign(taxonomy='224308', ncbi_taxonomy='224308')

        gene = create_input_value(df=gene, col='locus_tag')
        return gene

    def transform(self):
        gene = read(source=self.source, version=self.version, file='Gene.json', reader=read_json_lines,
                    default=pd.DataFrame(columns=['name', 'url', 'regulation', 'pubmed', 'tf', 'tfbs']))

        genome_path = Settings.genomes_database.joinpath(f'224308.json')
        genome = read_json_frame(genome_path)
        genome = genome[['locus_tag', 'name']].copy()

        genome = genome.assign(name_lower=genome['name'].str.lower())
        genome = genome.drop(columns=['name'])

        genes = self.transform_gene(gene=gene, genome=genome)

        ProtrendReporter.report_objects(source=self.source, version=self.version,
                                        system='extract', label=self.node.node_name(),
                                        objects=genes.shape[0], properties=genes.shape[1])

        annotated_genes = self.annotate_genes(genes)

        df = self.merge_annotations(annotated_genes, genes)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
