from typing import List, Union

import pandas as pd

from protrend.io.json import read_json_lines, read_json_frame
from protrend.io.utils import read_from_stack
from protrend.transform.annotation import annotate_genes
from protrend.transform.connector import DefaultConnector
from protrend.transform.dto import GeneDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, to_int_str
from protrend.transform.regprecise.organism import OrganismTransformer
from protrend.transform.regprecise.settings import RegulatorSettings, RegulatorToSource, RegulatorToOrganism
from protrend.transform.regprecise.source import SourceTransformer
from protrend.transform.transformer import Transformer


class RegulatorTransformer(Transformer):
    default_settings = RegulatorSettings
    columns = {'protrend_id',
               'organism_protrend_id', 'genome_id', 'ncbi_taxonomy',
               'mechanism',
               'regulon_id', 'name', 'genome', 'url', 'regulator_type', 'rfam',
               'biological_process', 'regulation_effector', 'regulation_regulog',
               'regulog', 'taxonomy', 'rna_family', 'effector', 'pathway', 'operon',
               'tfbs', 'gene', 'regulator_locus_tag', 'regulator_family',
               'regulation_mode', 'transcription_factor', 'tf_family',
               'locus_tag', 'synonyms', 'function', 'description',
               'ncbi_gene', 'ncbi_protein',
               'genbank_accession', 'refseq_accession', 'uniprot_accession',
               'sequence', 'strand', 'position_left', 'position_right', }

    read_columns = {'regulon_id', 'name', 'genome', 'url', 'regulator_type', 'rfam',
                    'biological_process', 'regulation_effector', 'regulation_regulog',
                    'regulog', 'taxonomy', 'rna_family', 'effector', 'pathway', 'operon',
                    'tfbs', 'gene', 'regulator_locus_tag', 'regulator_family',
                    'regulation_mode', 'transcription_factor', 'tf_family'}

    def _transform_tf(self, regulon: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:

        # filter tfs only
        regulon = regulon.dropna(subset=['regulator_locus_tag', 'genome'])

        regulon = self.drop_duplicates(df=regulon, subset=['regulator_locus_tag'], perfect_match=True,
                                       preserve_nan=False)

        apply_processors(rstrip, lstrip, df=regulon, col='regulator_locus_tag')
        apply_processors(rstrip, lstrip, df=regulon, col='name')

        regulon['mechanism'] = 'transcription factor'

        df = pd.merge(regulon, organism, left_on='genome', right_on='genome_id')

        self.create_input_value(df=df, col='regulator_locus_tag')
        return df

    def _transform_rna(self, regulon: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:

        # filter tfs only
        regulon = regulon.dropna(subset=['rfam'])

        regulon = self.drop_duplicates(df=regulon, subset=['rfam', 'genome'], perfect_match=True, preserve_nan=False)

        apply_processors(rstrip, lstrip, df=regulon, col='rfam')
        apply_processors(rstrip, lstrip, df=regulon, col='name')

        regulon['mechanism'] = 'small RNA (sRNA)'

        df = pd.merge(regulon, organism, left_on='genome', right_on='genome_id')

        self.create_input_value(df=df, col='name')
        return df

    @staticmethod
    def _annotate_genes(loci: List[Union[None, str]], names: List[str], taxa: List[int]):

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
        # position_left: List[int]
        # position_right: List[int]

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        regulon = read_from_stack(stack=self._transform_stack, file='regulon',
                                  default_columns=self.read_columns, reader=read_json_lines)

        apply_processors(to_int_str, df=regulon, col='regulon_id')
        apply_processors(to_int_str, df=regulon, col='genome')

        organism = read_from_stack(stack=self._transform_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = self.select_columns(organism, 'protrend_id', 'genome_id', 'ncbi_taxonomy')
        organism = organism.rename(columns={'protrend_id': 'organism_protrend_id'})
        apply_processors(to_int_str, df=organism, col='genome_id')
        apply_processors(to_int_str, df=organism, col='ncbi_taxonomy')

        # ------------------ regulon of type TF --------------------------------
        tf = self._transform_tf(regulon=regulon, organism=organism)

        loci = tf['input_value'].tolist()
        names = tf['name'].tolist()
        taxa = tf['ncbi_taxonomy'].tolist()

        tfs = self._annotate_genes(loci, names, taxa)

        tf = pd.merge(tfs, tf, on='input_value', suffixes=('_annotation', '_regprecise'))

        tf = self.merge_columns(df=tf, column='name', left='name_annotation', right='name_regprecise')

        tf['locus_tag'] = tf['locus_tag'].fillna(tf['regulator_locus_tag'])

        tf = tf.drop(columns=['input_value'])

        # ------------------ regulon of type RNA --------------------------------
        rna = self._transform_rna(regulon=regulon, organism=organism)

        loci = [None] * rna.shape[0]
        names = rna['name'].tolist()
        taxa = rna['ncbi_taxonomy'].tolist()

        rnas = self._annotate_genes(loci, names, taxa)

        rna = pd.merge(rnas, rna, on='input_value', suffixes=('_annotation', '_regprecise'))

        rna = self.merge_columns(df=rna, column='name', left='name_annotation', right='name_regprecise')

        rna = rna.drop(columns=['input_value'])

        # --------------------- concat DFs --------------------------------------
        df = pd.concat([tf, rna], axis=0)

        apply_processors(to_int_str, df=df, col='genome_id')
        apply_processors(to_int_str, df=df, col='ncbi_taxonomy')
        apply_processors(to_int_str, df=df, col='genome')
        apply_processors(to_int_str, df=df, col='regulog')

        self._stack_transformed_nodes(df)
        return df


class RegulatorToSourceConnector(DefaultConnector):
    default_settings = RegulatorToSource

    def connect(self):
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = regulator['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=regulator['url'].tolist(),
                      external_identifier=regulator['regulon_id'].tolist(),
                      key=['regulon_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_csv(df)


class RegulatorToOrganismConnector(DefaultConnector):
    default_settings = RegulatorToOrganism

    def connect(self):
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)

        from_identifiers = regulator['protrend_id'].tolist()
        to_identifiers = regulator['organism_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)
