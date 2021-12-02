from collections import defaultdict
from typing import List, Dict, Union

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from protrend.io import read_json_lines, read_from_stack
from protrend.model import Regulator
from protrend.annotation import annotate_genes, GeneDTO
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.utils.processors import rstrip, lstrip, apply_processors
from protrend.utils import SetList


class RegulatorTransformer(DBTBSTransformer,
                           source='dbtbs',
                           version='0.0.3',
                           node=Regulator,
                           order=100,
                           register=True):
    default_transform_stack = {'tf': 'TranscriptionFactor.json', 'sequence': 'sequence.gb'}
    columns = SetList(['locus_tag', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession',
                       'uniprot_accession', 'sequence', 'strand', 'start', 'stop', 'mechanism', 'protrend_id',
                       'name_dbtbs', 'family', 'domain', 'domain_description', 'url',
                       'type', 'comment', 'operon', 'subti_list', 'consensus_sequence'])

    read_columns = SetList(['name', 'family', 'domain', 'domain_description', 'description', 'url',
                            'type', 'comment', 'operon', 'subti_list', 'consensus_sequence'])

    @staticmethod
    def _get_gene_name(qualifiers: Dict) -> Union[str, None]:
        gene: List[str] = qualifiers.get('gene', [None])

        if gene[0]:
            return gene[0].lower()

        return

    @staticmethod
    def _get_locus(qualifiers: Dict) -> Union[str, None]:
        locus: List[str] = qualifiers.get('locus_tag', [None])

        if locus[0]:
            return locus[0].lower()

        return

    @staticmethod
    def _get_genbank(qualifiers: Dict) -> Union[str, None]:
        gb: List[str] = qualifiers.get('protein_id', [None])

        if gb[0]:
            return gb[0]

        return

    @staticmethod
    def _get_uniprot(qualifiers: Dict) -> Union[str, None]:
        xrefs: List[str] = qualifiers.get('db_xref', [])

        for xref in xrefs:

            if 'UniProtKB/Swiss-Prot:' in xref:
                return xref.replace('UniProtKB/Swiss-Prot:', '').rstrip().lstrip()

        return

    def _index_sequence(self, sequence: SeqRecord) -> pd.DataFrame:

        sequence_idx = defaultdict(SetList)

        for feature in sequence.features:

            if feature.type == 'CDS':

                gene_name = self._get_gene_name(feature.qualifiers)

                if gene_name:
                    sequence_idx[gene_name].append(gene_name)

                    locus = self._get_locus(feature.qualifiers)
                    sequence_idx[gene_name].append(locus)

                    gb_acc = self._get_genbank(feature.qualifiers)
                    sequence_idx[gene_name].append(gb_acc)

                    uniprot_acc = self._get_uniprot(feature.qualifiers)
                    sequence_idx[gene_name].append(uniprot_acc)

        # filter out duplicated gene names
        sequence_idx = {i: val for i, val in enumerate(sequence_idx.values())
                        if len(val) == 4}

        return pd.DataFrame.from_dict(sequence_idx,
                                      orient='index',
                                      columns=['gene_name_lower', 'locus_tag',
                                               'genbank_accession', 'uniprot_accession'])

    def _transform_tf(self, tf: pd.DataFrame, sequence: pd.DataFrame) -> pd.DataFrame:
        tf = tf.drop(columns=['description'])

        # filter nan and duplicates
        tf = self.drop_duplicates(df=tf, subset=['name'], perfect_match=True)
        tf = tf.dropna(subset=['name'])

        tf['mechanism'] = 'transcription factor'

        tf = apply_processors(tf, name=[rstrip, lstrip])

        tf['gene_name_lower'] = tf['name'].str.lower()

        tf = pd.merge(tf, sequence, on='gene_name_lower')

        tf = self.create_input_value(df=tf, col='locus_tag')
        return tf

    @staticmethod
    def _annotate_tfs(loci: List[str],
                      names: List[str],
                      genbanks: List[str],
                      accessions: List[str],
                      taxa: List[str]) -> pd.DataFrame:
        # annotate with uniprot_accessions, ncbi_proteins, locus_tag for ncbi_gene
        dtos = [GeneDTO(input_value=locus) for locus in loci]
        annotate_genes(dtos=dtos, loci=loci, names=names, taxa=taxa,
                       uniprot_proteins=accessions, ncbi_genbanks=genbanks)

        for dto, locus, name, acc, gb in zip(dtos, loci, names, accessions, genbanks):
            dto.synonyms.append(name)
            dto.synonyms.append(locus)

            dto.locus_tag.insert(0, locus)
            dto.name.insert(0, name)
            dto.uniprot_accession.insert(0, acc)
            dto.genbank_accession.insert(0, gb)

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

        tfs = pd.DataFrame([dto.to_dict() for dto in dtos])
        strand_mask = (tfs['strand'] != 'reverse') & (tfs['strand'] != 'forward')
        tfs.loc[strand_mask, 'strand'] = None
        return tfs

    def transform(self):
        tf = read_from_stack(stack=self.transform_stack, file='tf',
                             default_columns=self.read_columns, reader=read_json_lines)

        gb_file = self.transform_stack['sequence']
        sequence = SeqIO.read(gb_file, "genbank")
        sequence = self._index_sequence(sequence)

        tf = self._transform_tf(tf=tf, sequence=sequence)

        loci = tf['input_value'].tolist()
        names = tf['name'].tolist()
        gbs = tf['genbank_accession'].tolist()
        accessions = tf['uniprot_accession'].tolist()
        taxa = ['224308'] * len(loci)
        genes = self._annotate_tfs(loci=loci, names=names, genbanks=gbs, accessions=accessions, taxa=taxa)

        df = pd.merge(genes, tf, on='input_value', suffixes=('_annotation', '_dbtbs'))

        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_dbtbs')

        df['old_name'] = df['name_dbtbs']
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_dbtbs')
        df = df.rename(columns={'old_name': 'name_dbtbs'})

        df = self.merge_columns(df=df, column='genbank_accession',
                                left='genbank_accession_annotation', right='genbank_accession_dbtbs')
        df = self.merge_columns(df=df, column='uniprot_accession',
                                left='uniprot_accession_annotation', right='uniprot_accession_dbtbs')

        df = df.dropna(subset=['locus_tag'])

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)

        return df
