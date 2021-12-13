from collections import defaultdict
from typing import List, Dict, Union

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from protrend.io import read_json_lines, read_from_stack
from protrend.model import Gene
from protrend.annotation import annotate_genes, GeneDTO
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.utils.processors import (rstrip, lstrip, apply_processors)
from protrend.utils import SetList


class GeneTransformer(DBTBSTransformer,
                      source='dbtbs',
                      version='0.0.3',
                      node=Gene,
                      order=100,
                      register=True):
    default_transform_stack = {'gene': 'Gene.json', 'sequence': 'sequence.gb'}
    columns = SetList(['protrend_id', 'locus_tag', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession',
                       'uniprot_accession', 'sequence', 'strand', 'start', 'stop',
                       'url', 'position', 'cog_id', 'conversed_groups', 'tf',
                       'operon', 'tfbs', 'name_dbtbs'])
    read_columns = SetList(['name', 'url', 'synonyms', 'strand', 'position', 'cog_id', 'conversed_groups', 'tf',
                            'operon', 'tfbs', 'function'])

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

    def _transform_gene(self, gene: pd.DataFrame, sequence: pd.DataFrame) -> pd.DataFrame:
        gene = gene.drop(columns=['synonyms', 'strand', 'function'])

        gene = gene.explode(column='name')

        # filter nan and duplicates
        gene = self.drop_duplicates(df=gene, subset=['name'], perfect_match=True)
        gene = gene.dropna(subset=['name'])

        gene = apply_processors(gene, name=[rstrip, lstrip])

        gene['gene_name_lower'] = gene['name'].str.lower()

        gene = pd.merge(gene, sequence, on='gene_name_lower')

        gene = self.create_input_value(df=gene, col='locus_tag')
        return gene

    @staticmethod
    def _annotate_genes(loci: List[str],
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

        genes = pd.DataFrame([dto.to_dict() for dto in dtos])
        strand_mask = (genes['strand'] != 'reverse') & (genes['strand'] != 'forward')
        genes.loc[strand_mask, 'strand'] = None
        return genes

    def transform(self):
        gene = read_from_stack(stack=self.transform_stack, key='gene',
                               columns=self.read_columns, reader=read_json_lines)

        gb_file = self.transform_stack['sequence']
        sequence = SeqIO.read(gb_file, "genbank")
        sequence = self._index_sequence(sequence)

        gene = self._transform_gene(gene=gene, sequence=sequence)

        loci = gene['input_value'].tolist()
        names = gene['name'].tolist()
        gbs = gene['genbank_accession'].tolist()
        accessions = gene['uniprot_accession'].tolist()
        taxa = ['224308'] * len(loci)
        genes = self._annotate_genes(loci=loci, names=names, genbanks=gbs, accessions=accessions, taxa=taxa)

        df = pd.merge(genes, gene, on='input_value', suffixes=('_annotation', '_dbtbs'))

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
