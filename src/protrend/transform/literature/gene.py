import pandas as pd

from protrend.model import Gene
from protrend.transform.literature.base import LiteratureTransformer, read_literature_networks
from protrend.transform.mix_ins import GeneMixIn
from protrend.utils import SetList


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

        self.stack_transformed_nodes(df)
        return df
