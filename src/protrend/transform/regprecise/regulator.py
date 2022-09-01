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
from protrend.utils.processors import rstrip, lstrip, apply_processors, to_int_str


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

    def transform_regulator(self, regulon: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        regulon = regulon.dropna(subset=['genome'])
        regulon = drop_empty_string(regulon, 'genome')
        regulon = apply_processors(regulon, regulon_id=to_int_str, genome=to_int_str)

        # + "ncbi_taxonomy"
        regulon = pd.merge(regulon, organism, on='genome')

        ProtrendReporter.report_objects(source=self.source, version=self.version,
                                        system='extract', label=self.node.node_name(),
                                        objects=regulon.shape[0], properties=regulon.shape[1])

        # mechanism
        mechanism_map = {'RNA regulatory element': SMALL_RNA, 'Transcription factor': TRANSCRIPTION_FACTOR}
        mechanism = regulon['regulator_type'].map(mechanism_map)

        regulon = regulon.assign(locus_tag=regulon['regulator_locus_tag'].copy(),
                                 regulator_name=regulon['name'].copy(),
                                 mechanism=mechanism)
        regulon = apply_processors(regulon, locus_tag=[rstrip, lstrip], name=[rstrip, lstrip])

        # select transcription factors only
        regulon = regulon[regulon['mechanism'] == TRANSCRIPTION_FACTOR]
        regulon = regulon.reset_index(drop=True)

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

        annotated_regulators = self.annotate_genes(regulators)

        df = pd.merge(annotated_regulators, regulators, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_regprecise')
        df = merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        df = df.dropna(subset=['locus_tag'])
        df = drop_empty_string(df, 'locus_tag')
        df = drop_duplicates(df, subset=['locus_tag'])

        df = apply_processors(df, regulon_id=to_int_str, genome=to_int_str, ncbi_taxonomy=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
