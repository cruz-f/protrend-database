import pandas as pd

from protrend.io import read_json_lines, read_json_frame, read_from_stack
from protrend.model import Regulator
from protrend.transform.mix_ins import GeneMixIn
from protrend.transform.regprecise.base import RegPreciseTransformer
from protrend.transform.regprecise.organism import OrganismTransformer
from protrend.transform.transformations import drop_empty_string, create_input_value, merge_columns
from protrend.utils import SetList
from protrend.utils.processors import rstrip, lstrip, apply_processors, to_int_str


class RegulatorTransformer(GeneMixIn, RegPreciseTransformer,
                           source='regprecise',
                           version='0.0.0',
                           node=Regulator,
                           order=90,
                           register=True):
    default_transform_stack = {'organism': 'integrated_organism.json', 'regulon': 'Regulon.json'}
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop', 'mechanism',
                       'ncbi_taxonomy', 'regulator_name', 'regulator_locus_tag',
                       'regulon_id', 'genome', 'url', 'regulator_type', 'rfam',
                       'regulator_family', 'regulation_mode', 'biological_process', 'regulation_effector',
                       'regulation_regulog', 'regulog', 'taxonomy', 'transcription_factor', 'tf_family',
                       'rna_family', 'effector', 'pathway', 'operon', 'tfbs', 'gene'])
    read_columns = SetList(['regulon_id', 'name', 'genome', 'url', 'regulator_type', 'rfam', 'regulator_locus_tag',
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
        mechanism_map = {'RNA regulatory element': 'small RNA (sRNA)', 'Transcription factor': 'transcription factor'}
        mechanism = regulon['regulator_type'].map(mechanism_map)

        regulon = regulon.assign(locus_tag=regulon['regulator_locus_tag'].copy(),
                                 regulator_name=regulon['name'].copy(),
                                 mechanism=mechanism)
        regulon = apply_processors(regulon, locus_tag=[rstrip, lstrip], name=[rstrip, lstrip])

        regulon = create_input_value(df=regulon, col='regulon_id')
        return regulon

    def transform(self):
        regulon = read_from_stack(stack=self.transform_stack, key='regulon',
                                  columns=self.read_columns, reader=read_json_lines)

        organism = read_from_stack(stack=self.transform_stack, key='organism',
                                   columns=OrganismTransformer.columns, reader=read_json_frame)

        organism = self.transform_organism(organism)
        regulators = self.transform_regulator(regulon, organism)

        annotated_regulators = self.annotate_genes(regulators)

        df = pd.merge(annotated_regulators, regulators, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        df = merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_regprecise')

        # the small RNAs might not have any locus tag associated with during the annotation, so we will create new
        # locus tag composed by the name of sRNA plus the taxonomy identifier
        loci = df['name'] + df['ncbi_taxonomy']
        loci = df['locus_tag'].fillna(loci)
        df = df.assign(locus_tag=loci)

        df = apply_processors(df, regulon_id=to_int_str, genome=to_int_str, ncbi_taxonomy=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
