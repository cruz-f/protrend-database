import numpy as np
import pandas as pd

from protrend.io import read_json_lines, read
from protrend.io.utils import read_organism
from protrend.model import Regulator
from protrend.report import ProtrendReporter
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.regprecise.base import RegPreciseTransformer
from protrend.transform.regprecise.organism import OrganismTransformer
from protrend.transform.transformations import drop_empty_string, create_input_value, merge_columns, drop_duplicates
from protrend.utils import SetList
from protrend.utils.constants import SMALL_RNA, TRANSCRIPTION_FACTOR
from protrend.utils.processors import rstrip, lstrip, apply_processors, to_int_str, to_str


class RegulatorTransformer(GeneMixIn, RegPreciseTransformer,
                           source='regprecise',
                           version='0.0.0',
                           node=Regulator,
                           order=90,
                           register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'protein_sequence', 'strand', 'start', 'stop', 'mechanism',
                       'ncbi_taxonomy', 'regulator_name', 'regulator_locus_tag',
                       'regulon_id', 'genome', 'url', 'regulator_type', 'rfam',
                       'regulator_family', 'regulation_mode', 'biological_process', 'regulation_effector',
                       'regulation_regulog', 'regulog', 'taxonomy', 'transcription_factor', 'tf_family',
                       'rna_family', 'effector', 'pathway', 'operon', 'tfbs', 'gene'])

    @staticmethod
    def transform_regulator(regulon: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        regulon = regulon.dropna(subset=['genome'])
        regulon = drop_empty_string(regulon, 'genome')
        regulon = apply_processors(regulon, regulon_id=to_int_str, genome=to_int_str)

        # + "ncbi_taxonomy"
        regulon = pd.merge(regulon, organism, on='genome')

        # mechanism
        mechanism_map = {'RNA regulatory element': SMALL_RNA, 'Transcription factor': TRANSCRIPTION_FACTOR}
        mechanism = regulon['regulator_type'].map(mechanism_map)

        regulon = regulon.assign(locus_tag=regulon['regulator_locus_tag'].copy(),
                                 regulator_name=regulon['name'].copy(),
                                 mechanism=mechanism)
        regulon = apply_processors(regulon, locus_tag=[rstrip, lstrip], name=[rstrip, lstrip])

        regulon = create_input_value(df=regulon, col='regulon_id')
        return regulon

    def transform(self):
        regulon = read(source=self.source, version=self.version,
                       file='Regulon.json', reader=read_json_lines,
                       default=pd.DataFrame(columns=['regulon_id', 'name', 'genome', 'url', 'regulator_type', 'rfam',
                                                     'regulator_locus_tag',
                                                     'regulator_family', 'regulation_mode', 'biological_process',
                                                     'regulation_effector',
                                                     'regulation_regulog', 'regulog', 'taxonomy',
                                                     'transcription_factor', 'tf_family',
                                                     'rna_family', 'effector', 'pathway', 'operon', 'tfbs', 'gene']))

        organism = read_organism(source=self.source, version=self.version, columns=OrganismTransformer.columns)

        organism = self.transform_organism(organism)
        regulators = self.transform_regulator(regulon, organism)

        ProtrendReporter.report_objects(source=self.source, version=self.version,
                                        system='extract', label=self.node.node_name(),
                                        objects=regulators.shape[0], properties=regulators.shape[1])

        annotated_regulators = self.annotate_genes(regulators, drop_nan_locus_tag=False)

        df = pd.merge(annotated_regulators, regulators, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_regprecise')
        loci = df['locus_tag'].replace('', np.nan)
        df = df.assign(locus_tag=loci)

        df = merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        # the small RNAs might not have any locus tag associated with during the annotation, so we will create new
        # locus tag composed by the name of sRNA plus the taxonomy identifier
        fake_ncbi = df['ncbi_taxonomy'].copy()
        fake_name = df['name'].copy()
        fake_str = '_'
        df = df.assign(fake_name=fake_name, fake_str=fake_str, fake_ncbi=fake_ncbi)
        df = apply_processors(df, fake_name=to_str, fake_str=to_str, fake_ncbi=to_int_str)
        fake_loci = df['fake_name'] + df['fake_str'] + df['fake_ncbi']

        loci = df['locus_tag'].fillna(fake_loci)
        df = df.assign(locus_tag=loci)
        df = df.drop(columns=['fake_name', 'fake_str', 'fake_ncbi'])

        df = df.dropna(subset=['locus_tag'])
        df = drop_empty_string(df, 'locus_tag')
        df = drop_duplicates(df, subset=['locus_tag'])

        df = apply_processors(df, regulon_id=to_int_str, genome=to_int_str, ncbi_taxonomy=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
