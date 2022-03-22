import numpy as np
import pandas as pd

from protrend import GeneDTO, annotate_genes
from protrend.io import read_json_lines, read
from protrend.io.utils import read_organism
from protrend.log import ProtrendLogger
from protrend.model import Regulator
from protrend.transform.mix_ins import GeneMixIn, get_values
from protrend.transform.regprecise.base import RegPreciseTransformer
from protrend.transform.regprecise.organism import OrganismTransformer
from protrend.transform.transformations import drop_empty_string, create_input_value, merge_columns, drop_duplicates
from protrend.utils import SetList
from protrend.utils.constants import SMALL_RNA, TRANSCRIPTION_FACTOR, REVERSE, FORWARD, UNKNOWN
from protrend.utils.processors import rstrip, lstrip, apply_processors, to_int_str, to_str


class RegulatorTransformer(GeneMixIn, RegPreciseTransformer,
                           source='regprecise',
                           version='0.0.0',
                           node=Regulator,
                           order=90,
                           register=True):
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop', 'mechanism',
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

    @staticmethod
    def annotate_genes(df: pd.DataFrame) -> pd.DataFrame:
        input_values = get_values(df, 'input_value')

        genes = [GeneDTO(input_value=input_value) for input_value in input_values]

        ProtrendLogger.log.info(f'Annotating {len(genes)} genes')

        loci = get_values(df, 'locus_tag')
        names = get_values(df, 'name')
        taxa = get_values(df, 'ncbi_taxonomy')
        uniprot_proteins = get_values(df, 'uniprot_accession')
        ncbi_proteins = get_values(df, 'ncbi_protein')
        ncbi_genbanks = get_values(df, 'genbank_accession')
        ncbi_refseqs = get_values(df, 'refseq_accession')
        ncbi_genes = get_values(df, 'ncbi_gene')

        iterator = zip(
            ('locus_tag', 'name', 'ncbi_taxonomy', 'uniprot_accession', 'ncbi_protein', 'genbank_accession',
             'refseq_accession', 'ncbi_gene'),
            (loci, names, taxa, uniprot_proteins, ncbi_proteins, ncbi_genbanks, ncbi_refseqs, ncbi_genes)
        )

        params = [param for param, value in iterator if value is not None]
        params = ','.join(params)

        ProtrendLogger.log.info(f'Annotating with the following params: {params}')

        annotate_genes(dtos=genes,
                       loci=loci,
                       names=names,
                       taxa=taxa,
                       uniprot_proteins=uniprot_proteins,
                       ncbi_proteins=ncbi_proteins,
                       ncbi_genbanks=ncbi_genbanks,
                       ncbi_refseqs=ncbi_refseqs,
                       ncbi_genes=ncbi_genes)

        genes_dict = [dto.to_dict() for dto in genes]
        genes_df = pd.DataFrame(genes_dict)

        if genes_df.empty:
            return genes_df

        strand_mask = (genes_df['strand'] != REVERSE) & (genes_df['strand'] != FORWARD)
        genes_df.loc[strand_mask, 'strand'] = UNKNOWN

        return genes_df

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
