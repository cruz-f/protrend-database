from typing import List, Union

import pandas as pd

from protrend.io import read_json_lines, read_json_frame, read_from_stack
from protrend.model import Regulator, Source, Organism
from protrend.annotation import annotate_genes, GeneDTO
from protrend.utils.processors import rstrip, lstrip, apply_processors, to_int_str
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.regprecise.organism import OrganismTransformer
from protrend.transform.regprecise.source import SourceTransformer
from protrend.utils import SetList


class RegulatorTransformer(RegPreciseTransformer,
                           source='regprecise',
                           version='0.0.0',
                           node=Regulator,
                           order=90,
                           register=True):
    default_transform_stack = {'regulon': 'Regulon.json',
                               'organism': 'integrated_organism.json'}
    columns = SetList(['locus_tag', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession',
                       'uniprot_accession', 'sequence', 'strand', 'start', 'stop',
                       'regulon_id', 'genome', 'url', 'regulator_type', 'rfam',
                       'biological_process', 'regulation_effector', 'regulation_regulog',
                       'regulog', 'taxonomy', 'rna_family', 'effector', 'pathway', 'operon',
                       'tfbs', 'gene', 'regulator_locus_tag', 'regulator_family',
                       'regulation_mode', 'transcription_factor', 'tf_family', 'mechanism',
                       'organism_protrend_id', 'genome_id', 'ncbi_taxonomy', 'name',
                       'protrend_id'])

    read_columns = SetList(['regulon_id', 'name', 'genome', 'url', 'regulator_type', 'rfam', 'biological_process',
                            'regulation_effector', 'regulation_regulog', 'regulog', 'taxonomy', 'rna_family',
                            'effector', 'pathway', 'operon', 'tfbs', 'gene'])

    def _transform_tf(self, regulon: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        # filter tfs only
        regulon = regulon.dropna(subset=['regulator_locus_tag', 'genome'])

        regulon = self.drop_duplicates(df=regulon, subset=['regulator_locus_tag'], perfect_match=True,
                                       preserve_nan=False)

        regulon = apply_processors(regulon, regulator_locus_tag=[rstrip, lstrip], name=[rstrip, lstrip])

        regulon['mechanism'] = 'transcription factor'

        df = pd.merge(regulon, organism, how='left', left_on='genome', right_on='genome_id')

        self.create_input_value(df=df, col='regulator_locus_tag')
        return df

    def _transform_rna(self, regulon: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        # filter tfs only
        regulon = regulon.dropna(subset=['rfam', 'genome'])

        regulon = self.drop_duplicates(df=regulon, subset=['rfam', 'genome'], perfect_match=True, preserve_nan=False)

        regulon = apply_processors(regulon, rfam=[rstrip, lstrip], name=[rstrip, lstrip])

        regulon['mechanism'] = 'small RNA (sRNA)'

        df = pd.merge(regulon, organism, how='left', left_on='genome', right_on='genome_id')

        return df

    @staticmethod
    def _annotate_tfs(loci: List[Union[None, str]], names: List[str], taxa: List[str]):
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

        regulators = pd.DataFrame([dto.to_dict() for dto in dtos])
        strand_mask = (regulators['strand'] != 'reverse') & (regulators['strand'] != 'forward')
        regulators.loc[strand_mask, 'strand'] = None
        return regulators

    @staticmethod
    def _annotate_rnas(names: List[str]):
        dtos = [GeneDTO(input_value=name) for name in names]

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

        df = pd.DataFrame([dto.to_dict() for dto in dtos])
        df = df.drop(columns=['name', 'locus_tag', 'input_value'])
        return df

    def transform(self):
        regulon = read_from_stack(stack=self.transform_stack, key='regulon',
                                  columns=self.read_columns, reader=read_json_lines)

        regulon = apply_processors(regulon, regulon_id=to_int_str, genome=to_int_str)

        organism = read_from_stack(stack=self.transform_stack, key='organism',
                                   columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = self.select_columns(organism, 'protrend_id', 'genome_id', 'ncbi_taxonomy')
        organism = organism.rename(columns={'protrend_id': 'organism_protrend_id'})

        organism = apply_processors(organism, genome_id=to_int_str, ncbi_taxonomy=to_int_str)

        # ------------------ regulon of type TF --------------------------------
        tf = self._transform_tf(regulon=regulon, organism=organism)

        loci = tf['input_value'].tolist()
        names = tf['name'].tolist()
        taxa = tf['ncbi_taxonomy'].tolist()

        tfs = self._annotate_tfs(loci, names, taxa)

        tf = pd.merge(tfs, tf, on='input_value', suffixes=('_annotation', '_regprecise'))

        tf = self.merge_columns(df=tf, column='name', left='name_annotation', right='name_regprecise')

        tf['locus_tag'] = tf['locus_tag'].fillna(tf['regulator_locus_tag'])

        tf = tf.drop(columns=['input_value'])

        # ------------------ regulon of type RNA --------------------------------
        rna = self._transform_rna(regulon=regulon, organism=organism)

        names = rna['name'].tolist()
        rnas = self._annotate_rnas(names)

        rna = pd.concat([rna, rnas], axis=1)
        rna['locus_tag'] = rna['organism_protrend_id'] + '_' + rna['rfam']

        # --------------------- concat DFs --------------------------------------
        df = pd.concat([tf, rna])

        df = apply_processors(df, genome_id=to_int_str, ncbi_taxonomy=to_int_str, genome=to_int_str, regulog=to_int_str)

        self.stack_transformed_nodes(df)
        return df


class RegulatorToSourceConnector(RegPreciseConnector,
                                 source='regprecise',
                                 version='0.0.0',
                                 from_node=Regulator,
                                 to_node=Source,
                                 register=True):
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'source': 'integrated_source.json'}

    def connect(self):
        regulator = read_from_stack(stack=self._connect_stack, key='regulator',
                                    columns=RegulatorTransformer.columns, reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, key='source',
                                 columns=SourceTransformer.columns, reader=read_json_frame)

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

        self.stack_json(df)


class RegulatorToOrganismConnector(RegPreciseConnector,
                                   source='regprecise',
                                   version='0.0.0',
                                   from_node=Regulator,
                                   to_node=Organism,
                                   register=True):
    default_connect_stack = {'regulator': 'integrated_regulator.json'}

    def connect(self):
        regulator = read_from_stack(stack=self._connect_stack, key='regulator',
                                    columns=RegulatorTransformer.columns, reader=read_json_frame)

        from_identifiers = regulator['protrend_id'].tolist()
        to_identifiers = regulator['organism_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
