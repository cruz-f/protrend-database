import pandas as pd

from protrend.io import read_from_multi_stack
from protrend.model import Gene
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip


class GeneTransformer(CoryneRegNetTransformer,
                      source='coryneregnet',
                      version='0.0.0',
                      node=Gene,
                      order=100,
                      register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop',
                       'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                       'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                       'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence',
                       'PMID', 'Source', 'taxonomy', 'source'])

    def transform_gene(self, network: pd.DataFrame) -> pd.DataFrame:
        gene = network.dropna(subset=['TG_locusTag'])
        gene = self.drop_empty_string(gene, 'TG_locusTag')
        gene = self.drop_duplicates(df=gene, subset=['TG_locusTag'])

        gene = apply_processors(gene, TG_locusTag=[rstrip, lstrip], TG_name=[rstrip, lstrip])

        gene = gene.assign(locus_tag=gene['TG_locusTag'].copy(), name=gene['TG_name'].copy(),
                           ncbi_taxonomy=gene['taxonomy'].copy())

        gene = self.create_input_value(df=gene, col='locus_tag')
        return gene

    def transform(self):
        network = read_from_multi_stack(stack=self.transform_stack, key='network', columns=self.default_network_columns)

        genes = self.transform_gene(network)
        annotated_genes = self.annotate_genes(genes)

        df = self.merge_annotations(annotated_genes, genes)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
