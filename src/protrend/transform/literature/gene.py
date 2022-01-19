import pandas as pd

from protrend.model import Gene
from protrend.transform.literature.base import LiteratureTransformer, read_literature_networks
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.transformations import drop_empty_string, drop_duplicates
from protrend.utils import SetList, apply_processors
from protrend.utils.processors import to_str, to_int_str


class GeneTransformer(GeneMixIn, LiteratureTransformer,
                      source='literature',
                      version='0.0.0',
                      node=Gene,
                      order=100,
                      register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop',
                       'regulator_locus_tag', 'gene_locus_tag',
                       'regulatory_effect', 'evidence', 'effector_name', 'mechanism',
                       'publication', 'taxonomy', 'source'])

    def transform_gene(self, network: pd.DataFrame) -> pd.DataFrame:
        return self._transform_gene(network, col='gene_locus_tag')

    def transform(self):
        network = read_literature_networks(source=self.source, version=self.version)

        genes = self.transform_gene(network)
        annotated_genes = self.annotate_genes(genes)

        df = self.merge_annotations(annotated_genes, genes)

        # the small RNAs might not have any locus tag associated with during the annotation, so we will create new
        # locus tag composed by the name of sRNA plus the taxonomy identifier
        fake_ncbi = df['taxonomy'].copy()
        fake_name = df['name'].copy()
        fake_str = ' for organism '
        df = df.assign(fake_name=fake_name, fake_str=fake_str, fake_ncbi=fake_ncbi)
        df = apply_processors(df, fake_name=to_str, fake_str=to_str, fake_ncbi=to_int_str)
        fake_loci = df['fake_name'] + df['fake_str'] + df['fake_ncbi']

        loci = df['locus_tag'].fillna(fake_loci)
        df = df.assign(locus_tag=loci)
        df = df.drop(columns=['fake_name', 'fake_str', 'fake_ncbi'])

        df = df.dropna(subset=['locus_tag'])
        df = drop_empty_string(df, 'locus_tag')
        df = drop_duplicates(df, subset=['locus_tag'])

        self.stack_transformed_nodes(df)
        return df
