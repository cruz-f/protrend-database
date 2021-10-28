from typing import List, Union, Tuple

import pandas as pd

from protrend.model.model import Regulator
from protrend.transform import GeneDTO
from protrend.transform.annotation import annotate_genes
from protrend.transform.literature.base import LiteratureTransformer
from protrend.transform.processors import apply_processors, rstrip, lstrip, to_set_list
from protrend.utils import SetList


class GeneTransformer(LiteratureTransformer):
    default_node = Regulator
    default_order = 100
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession', 'sequence',
                       'strand', 'start', 'stop', 'mechanism',
                       'regulator_locus_tag', 'regulator_name', 'operon', 'genes_locus_tag',
                       'genes_name', 'regulatory_effect', 'evidence', 'effector', 'mechanism',
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

    def _transform_regulator(self, network: pd.DataFrame) -> pd.DataFrame:
        network = apply_processors(network, regulator_locus_tag=[rstrip, lstrip], regulator_name=[rstrip, lstrip])

        network = self.drop_duplicates(df=network, subset=['regulator_locus_tag'],
                                       perfect_match=True, preserve_nan=True)
        network['regulator_locus_tag'] = network['regulator_locus_tag'].fillna(network['regulator_name'])
        network = network.dropna(subset=['regulator_locus_tag'])

        network['locus_tag'] = network['regulator_locus_tag']
        network['name'] = network['regulator_name']

        regulator_mechanisms = {key.rstrip().lstrip().lower(): value
                                for key, value in self.regulator_mechanisms.items()}

        def map_filter_mechanism(item: str) -> str:

            if item.rstrip().lstrip().lower() in regulator_mechanisms:
                return regulator_mechanisms[item]

            return 'unknown'

        network = apply_processors(network, mechanism=map_filter_mechanism)

        network = self.create_input_value(df=network, col='locus_tag')
        return network

    @staticmethod
    def _annotate_regulators(loci: List[Union[None, str]], names: List[str], taxa: List[str]):
        dtos = [GeneDTO(input_value=locus) for locus in loci]
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

        loci = regulator['input_value'].tolist()
        names = regulator['genes_name'].tolist()
        taxa = regulator['taxonomy'].tolist()

        regulators = self._annotate_regulators(loci, names, taxa)

        df = pd.merge(regulators, regulator, on='input_value', suffixes=('_annotation', '_literature'))

        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_literature')
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_literature')

        df = df.drop(columns=['input_value'])

        self._stack_transformed_nodes(df)

        return df
