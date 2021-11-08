from typing import List, Union

import pandas as pd

from protrend.model import Gene
from protrend.annotation import annotate_genes, GeneDTO
from protrend.transform.literature.base import LiteratureTransformer
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_set_list
from protrend.utils import SetList


class GeneTransformer(LiteratureTransformer,
                      source='literature',
                      version='0.0.0',
                      node=Gene,
                      order=100,
                      register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession', 'sequence',
                       'strand', 'start', 'stop',
                       'regulator_locus_tag', 'operon', 'genes_locus_tag',
                       'regulatory_effect', 'evidence', 'effector', 'mechanism',
                       'publication', 'taxonomy', 'source', 'network_id'])

    @staticmethod
    def _filter_ecol_locus_genes(df: pd.DataFrame) -> pd.Series:
        df = df.copy()
        mask = (df['genes_locus_tag'].str.startswith('b')) & (df['source'] == 'ecol_fang_et_al_2017')
        df = df.loc[mask, :]
        df['locus_tag'] = df['genes_locus_tag']
        df['name'] = None
        return df

    @staticmethod
    def _filter_mtub_locus_genes(df: pd.DataFrame) -> pd.Series:
        mask = (df['genes_locus_tag'].str.startswith('R')) & (df['source'] == 'mtub_turkarslan_et_al_2015')
        df = df.loc[mask, :]
        df['locus_tag'] = df['genes_locus_tag']
        df['name'] = None
        return df

    @staticmethod
    def _filter_paer_locus_genes(df: pd.DataFrame) -> pd.Series:
        mask = (df['genes_locus_tag'].str.startswith('PA')) & (df['source'] == 'paer_vasquez_et_al_2011')
        df = df.loc[mask, :]
        df['locus_tag'] = df['genes_locus_tag']
        df['name'] = None
        return df

    @staticmethod
    def _filter_paer_names_genes(df: pd.DataFrame) -> pd.Series:
        mask = (~ df['genes_locus_tag'].str.startswith('PA')) & (df['source'] == 'paer_vasquez_et_al_2011')
        df = df.loc[mask, :]
        df['locus_tag'] = None
        df['name'] = df['genes_locus_tag']
        return df

    @staticmethod
    def _filter_bsub_locus_genes(df: pd.DataFrame) -> pd.Series:
        mask = (df['genes_locus_tag'].str.startswith('BSU')) & (df['source'] == 'bsub_faria_et_al_2017')
        df = df.loc[mask, :]
        df['locus_tag'] = df['genes_locus_tag']
        df['name'] = None
        return df

    def _transform_gene(self, network: pd.DataFrame) -> pd.DataFrame:
        network = apply_processors(network, genes_locus_tag=to_set_list)
        network = network.explode(column='genes_locus_tag')

        network = apply_processors(network, genes_locus_tag=[rstrip, lstrip])
        network = self.drop_duplicates(df=network, subset=['genes_locus_tag', 'taxonomy'], perfect_match=True)
        network = network.dropna(subset=['genes_locus_tag'])

        filtered_networks = [self._filter_ecol_locus_genes(network),
                             self._filter_mtub_locus_genes(network),
                             self._filter_paer_locus_genes(network),
                             self._filter_paer_names_genes(network),
                             self._filter_bsub_locus_genes(network)]
        network = pd.concat(filtered_networks)
        network = network.reset_index(drop=True)

        network['gene_network_id'] = network['genes_locus_tag'] + network['taxonomy']

        network = self.create_input_value(df=network, col='gene_network_id')
        return network

    @staticmethod
    def _annotate_genes(input_values: Union[List[str], List[None]],
                        loci: List[Union[None, str]],
                        names: List[str],
                        taxa: List[str]):
        dtos = [GeneDTO(input_value=input_value) for input_value in input_values]
        annotate_genes(dtos=dtos, loci=loci, names=names, taxa=taxa)

        for dto, name in zip(dtos, names):
            dto.synonyms.append(name)

        # locus_tag: List[str]
        # name: List[str]
        # synonyms: List[str]
        # function: List[str]
        # description: List[str]
        # ncbi_gene: List[str]
        # ncbi_protein: List[str]
        # genbank_accession: List[str]
        # refseq_accession: List[str]
        # uniprot_accession: List[str]
        # sequence: List[str]
        # strand: List[str]
        # start: List[int]
        # stop: List[int]

        genes = pd.DataFrame([dto.to_dict() for dto in dtos])
        strand_mask = (genes['strand'] != 'reverse') & (genes['strand'] != 'forward')
        genes.loc[strand_mask, 'strand'] = None
        return genes

    def transform(self):
        network = self._build_network()
        gene = self._transform_gene(network)

        input_values = gene['input_value'].tolist()
        loci = gene['locus_tag'].tolist()
        names = gene['name'].tolist()
        taxa = gene['taxonomy'].tolist()

        genes = self._annotate_genes(input_values, loci, names, taxa)

        df = pd.merge(genes, gene, on='input_value', suffixes=('_annotation', '_literature'))

        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_literature')
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_literature')

        df = df.drop(columns=['input_value', 'gene_network_id'])

        self._stack_transformed_nodes(df)

        return df
