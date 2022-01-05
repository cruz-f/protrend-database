import pandas as pd

from protrend.model import Regulator
from protrend.transform.literature.base import LiteratureTransformer
from protrend.utils import SetList
from protrend.utils.processors import apply_processors


class RegulatorTransformer(LiteratureTransformer,
                           source='literature',
                           version='0.0.0',
                           node=Regulator,
                           order=100,
                           register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop', 'mechanism',
                       'regulator_locus_tag', 'gene_locus_tag',
                       'regulatory_effect', 'evidence', 'effector_name',
                       'publication', 'taxonomy', 'source'])

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

    def transform_regulator(self, network: pd.DataFrame) -> pd.DataFrame:
        network = self._transform_gene(network, col='regulator_locus_tag')

        regulator_mechanisms = {key.rstrip().lstrip().lower(): value
                                for key, value in self.regulator_mechanisms.items()}

        def map_filter_mechanism(item: str) -> str:
            item = item.rstrip().lstrip().lower()

            if item in regulator_mechanisms:
                return regulator_mechanisms[item]

            return 'unknown'

        network = apply_processors(network, mechanism=map_filter_mechanism)
        return network

    def transform(self):
        network = self.read_network()

        regulators = self.transform_regulator(network)
        annotated_regulators = self.annotate_genes(regulators)

        df = self.merge_annotations(annotated_regulators, regulators)

        self.stack_transformed_nodes(df)
        return df
