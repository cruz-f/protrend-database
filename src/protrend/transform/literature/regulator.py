import pandas as pd

from protrend.model import Regulator
from protrend.transform.literature.base import LiteratureTransformer, read_literature_networks
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.transformations import drop_empty_string, drop_duplicates
from protrend.utils import SetList
from protrend.utils.constants import SIGMA_FACTOR, TRANSCRIPTION_FACTOR, UNKNOWN, TRANSCRIPTION_TERMINATOR, SMALL_RNA
from protrend.utils.processors import apply_processors, to_str, to_int_str


class RegulatorTransformer(GeneMixIn, LiteratureTransformer,
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

    regulator_mechanisms = {"sigma factor": SIGMA_FACTOR,
                            "sigma factor - ECF type": SIGMA_FACTOR,
                            "TF": TRANSCRIPTION_FACTOR,
                            "TF-TC": TRANSCRIPTION_FACTOR,
                            "TF+M": TRANSCRIPTION_FACTOR,
                            "TF+unk": TRANSCRIPTION_FACTOR,
                            "TF+PP": TRANSCRIPTION_FACTOR,
                            "TF+PP+M": TRANSCRIPTION_FACTOR,
                            "TF+P": TRANSCRIPTION_FACTOR,
                            "TF+P+M": TRANSCRIPTION_FACTOR,
                            "TF+S": TRANSCRIPTION_FACTOR,
                            "TF-proteolyse": TRANSCRIPTION_FACTOR,
                            "TF+unstability": TRANSCRIPTION_FACTOR,
                            "TF+M+unk": TRANSCRIPTION_FACTOR,
                            "TF+ADN-M": TRANSCRIPTION_FACTOR,
                            "P-PTC": UNKNOWN,
                            "P-AT+PTS ": TRANSCRIPTION_TERMINATOR,
                            "P-AT+M": TRANSCRIPTION_TERMINATOR,
                            "P": UNKNOWN,
                            "Riboswitch": UNKNOWN,
                            "RNA switch": UNKNOWN,
                            "Rna-BPA": UNKNOWN,
                            "small regulatory RNA": SMALL_RNA,
                            "small RNA - translation": SMALL_RNA,
                            "RNA-anti antiterminator": TRANSCRIPTION_TERMINATOR,
                            "anti-sense RNA": UNKNOWN,
                            "unk": UNKNOWN,
                            "Translation regulation": UNKNOWN,
                            "silico-TF+M": TRANSCRIPTION_FACTOR,
                            "silico-TF+unk": TRANSCRIPTION_FACTOR,
                            "silico RNA switch": SMALL_RNA,
                            "silico-riboswitch": SMALL_RNA,
                            "silico-TF+TC": TRANSCRIPTION_FACTOR}

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
        network = read_literature_networks(source=self.source, version=self.version)

        regulators = self.transform_regulator(network)
        annotated_regulators = self.annotate_genes(regulators)

        df = self.merge_annotations(annotated_regulators, regulators)

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
