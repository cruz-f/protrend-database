from typing import List, Union

import pandas as pd

from protrend.model import Regulator
from protrend.annotation import annotate_genes, GeneDTO
from protrend.transform.literature.base import LiteratureTransformer
from protrend.utils.processors import apply_processors, rstrip, lstrip
from protrend.utils import SetList


class RegulatorTransformer(LiteratureTransformer,
                           source='literature',
                           version='0.0.0',
                           node=Regulator,
                           order=100,
                           register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession', 'sequence',
                       'strand', 'start', 'stop', 'mechanism',
                       'regulator_locus_tag', 'operon', 'genes_locus_tag',
                       'regulatory_effect', 'evidence', 'effector', 'mechanism',
                       'publication', 'taxonomy', 'source', 'network_id'])

    regulator_mechanisms = {"sigma factor": "sigma factor",
                            "sigma factor - ECF type": "sigma factor",
                            "TF": "transcription factor",
                            "TF-TC": "transcription factor",
                            "TF+M": "transcription factor",
                            "TF+unk": "transcription factor",
                            "TF+PP": "transcription factor",
                            "TF+PP+M": "transcription factor",
                            "TF+P": "transcription factor",
                            "TF+P+M": "transcription factor",
                            "TF+S": "transcription factor",
                            "TF-proteolyse": "transcription factor",
                            "TF+unstability": "transcription factor",
                            "TF+M+unk": "transcription factor",
                            "TF+ADN-M": "transcription factor",
                            "P-PTC": "unknown",
                            "P-AT+PTS ": "transcription terminator",
                            "P-AT+M": "transcription terminator",
                            "P": "unknown",
                            "Riboswitch": "unknown",
                            "RNA switch": "unknown",
                            "Rna-BPA": "unknown",
                            "small regulatory RNA": "small RNA (sRNA)",
                            "small RNA - translation": "small RNA (sRNA)",
                            "RNA-anti antiterminator": "transcription terminator",
                            "anti-sense RNA": "unknown",
                            "unk": "unknown",
                            "Translation regulation": "unknown",
                            "silico-TF+M": "transcription factor",
                            "silico-TF+unk": "transcription factor",
                            "silico RNA switch": "small RNA (sRNA)",
                            "silico-riboswitch": "small RNA (sRNA)",
                            "silico-TF+TC": "transcription factor"}

    @staticmethod
    def _filter_ecol_locus_genes(df: pd.DataFrame) -> pd.Series:
        df = df.copy()
        mask = (df['regulator_locus_tag'].str.startswith('b')) & (df['source'] == 'ecol_fang_et_al_2017')
        df = df.loc[mask, :]
        df['locus_tag'] = df['regulator_locus_tag']
        df['name'] = None
        return df

    @staticmethod
    def _filter_mtub_locus_genes(df: pd.DataFrame) -> pd.Series:
        mask = (df['regulator_locus_tag'].str.startswith('R')) & (df['source'] == 'mtub_turkarslan_et_al_2015')
        df = df.loc[mask, :]
        df['locus_tag'] = df['regulator_locus_tag']
        df['name'] = None
        return df

    @staticmethod
    def _filter_paer_locus_genes(df: pd.DataFrame) -> pd.Series:
        mask = (df['regulator_locus_tag'].str.startswith('PA')) & (df['source'] == 'paer_vasquez_et_al_2011')
        df = df.loc[mask, :]
        df['locus_tag'] = df['regulator_locus_tag']
        df['name'] = None
        return df

    @staticmethod
    def _filter_paer_names_genes(df: pd.DataFrame) -> pd.Series:
        mask = (~ df['regulator_locus_tag'].str.startswith('PA')) & (df['source'] == 'paer_vasquez_et_al_2011')
        df = df.loc[mask, :]
        df['locus_tag'] = None
        df['name'] = df['regulator_locus_tag']
        return df

    @staticmethod
    def _filter_bsub_locus_genes(df: pd.DataFrame) -> pd.Series:
        mask = (df['regulator_locus_tag'].str.startswith('BSU')) & (df['source'] == 'bsub_faria_et_al_2017')
        df = df.loc[mask, :]
        df['locus_tag'] = df['regulator_locus_tag']
        df['name'] = None
        return df

    def _transform_regulator(self, network: pd.DataFrame) -> pd.DataFrame:
        network = apply_processors(network, regulator_locus_tag=[rstrip, lstrip])

        network = self.drop_duplicates(df=network, subset=['regulator_locus_tag', 'taxonomy'],
                                       perfect_match=True, preserve_nan=True)
        network = network.dropna(subset=['regulator_locus_tag'])

        filtered_networks = [self._filter_ecol_locus_genes(network),
                             self._filter_mtub_locus_genes(network),
                             self._filter_paer_locus_genes(network),
                             self._filter_paer_names_genes(network),
                             self._filter_bsub_locus_genes(network)]
        network = pd.concat(filtered_networks, axis=0)
        network = network.reset_index(drop=True)

        regulator_mechanisms = {key.rstrip().lstrip().lower(): value
                                for key, value in self.regulator_mechanisms.items()}

        def map_filter_mechanism(item: str) -> str:
            item = item.rstrip().lstrip().lower()

            if item in regulator_mechanisms:
                return regulator_mechanisms[item]

            return 'unknown'

        network = apply_processors(network, mechanism=map_filter_mechanism)

        network['regulator_network_id'] = network['regulator_locus_tag'] + network['taxonomy']

        network = self.create_input_value(df=network, col='regulator_network_id')
        return network

    @staticmethod
    def _annotate_regulators(input_values: Union[List[str], List[None]],
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
        regulator = self._transform_regulator(network)

        input_values = regulator['input_value'].tolist()
        loci = regulator['locus_tag'].tolist()
        names = regulator['name'].tolist()
        taxa = regulator['taxonomy'].tolist()

        regulators = self._annotate_regulators(input_values, loci, names, taxa)

        df = pd.merge(regulators, regulator, on='input_value', suffixes=('_annotation', '_literature'))

        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_literature')
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_literature')

        df = df.drop(columns=['input_value', 'regulator_network_id'])

        self._stack_transformed_nodes(df)

        return df
