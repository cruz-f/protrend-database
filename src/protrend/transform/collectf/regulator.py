from typing import List, Union

import pandas as pd
from Bio.SeqRecord import SeqRecord

from protrend.bioapis import map_uniprot_identifiers, fetch_uniprot_record
from protrend.io import read_from_stack, read_json_lines, read_json_frame
from protrend.model import Regulator
from protrend.transform.collectf.base import CollectfTransformer
from protrend.transform.collectf.organism import OrganismTransformer
from protrend.utils import SetList
from protrend.utils.processors import take_first, flatten_set_list, apply_processors, to_list_nan, rstrip, lstrip


def map_accession(acc: str, mapping: pd.DataFrame) -> Union[str, None]:
    mask = mapping['From'] == acc
    to = mapping.loc[mask, 'To'].to_list()

    if to:
        return to[0]

    return None


def uniprot_record_locus_tag(record: SeqRecord) -> Union[str, None]:
    annotations = getattr(record, 'annotations', {})

    if annotations:
        loci = annotations.get('gene_name_ordered locus', [])

        if not loci:
            loci = annotations.get('gene_name_ORF', [])

        for locus in loci:
            return locus

    return None


class RegulatorTransformer(CollectfTransformer,
                           source='collectf',
                           version='0.0.1',
                           node=Regulator,
                           order=90,
                           register=True):
    default_transform_stack = {'regulon': 'Regulon.json', 'organism': 'integrated_organism.json'}
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop', 'mechanism',
                       'url', 'organism', 'operon', 'gene', 'tfbs', 'experimental_evidence'
                       'organism_protrend_id', 'organism_name_collectf', 'ncbi_taxonomy'])
    read_columns = SetList(['uniprot_accession', 'name', 'url', 'organism', 'operon',
                            'gene', 'tfbs', 'experimental_evidence'])

    @staticmethod
    def get_ncbi_proteins_from_uniprot(uniprot_accessions: List[str]) -> List[Union[str, None]]:
        # map uniprot_accessions to ncbi_proteins
        uniprot_ncbi_proteins = map_uniprot_identifiers(uniprot_accessions, from_='ACC', to='P_GI')
        return [map_accession(accession, uniprot_ncbi_proteins) for accession in uniprot_accessions]

    @staticmethod
    def get_locus_tag_from_uniprot(uniprot_accessions: List[str]) -> List[Union[str, None]]:
        # fetch uniprot record to retrieve locus_tag
        uniprot_records = [fetch_uniprot_record(accession) for accession in uniprot_accessions]
        return [uniprot_record_locus_tag(record) for record in uniprot_records]

    def transform_regulon(self, regulon: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        regulon = apply_processors(regulon, uniprot_accession=[rstrip, lstrip], name=[rstrip, lstrip])
        regulon = regulon.dropna(subset=['uniprot_accession'])
        regulon = self.drop_empty_string(regulon, 'uniprot_accession')
        regulon = self.drop_duplicates(df=regulon, subset=['uniprot_accession'])

        df = pd.merge(regulon, organism, how='left', left_on='organism', right_on='organism_name_collectf')

        df = self.drop_duplicates(df=df, subset=['uniprot_accession', 'organism'], perfect_match=True)

        uniprot_accessions = df['uniprot_accession'].to_list()
        ncbi_proteins = self.get_ncbi_proteins_from_uniprot(uniprot_accessions)
        loci = self.get_locus_tag_from_uniprot(uniprot_accessions)

        df = df.assign(mechanism='transcription factor', ncbi_protein=ncbi_proteins, locus_tag=loci)

        df = self.create_input_value(df=df, col='uniprot_accession')
        return df

    def transform_organism(self, organism: pd.DataFrame):
        organism = self.select_columns(organism, 'protrend_id', 'name_collectf', 'ncbi_taxonomy')
        organism = organism.rename(columns={'protrend_id': 'organism_protrend_id',
                                            'name_collectf': 'organism_name_collectf'})
        return organism

    def transform(self):
        regulon = read_from_stack(stack=self.transform_stack, key='regulon',
                                  columns=self.read_columns, reader=read_json_lines)

        organism = read_from_stack(stack=self.transform_stack, key='organism',
                                   columns=OrganismTransformer.columns, reader=read_json_frame)

        organism = self.transform_organism(organism)
        regulators = self.transform_regulon(regulon, organism)
        annotated_regulators = self.annotate_genes(regulators)

        df = pd.merge(annotated_regulators, regulators, on='input_value', suffixes=('_annotation', '_collectf'))

        # merge loci
        df = self.merge_loci(df=df, left_suffix='_annotation', right_suffix='_collectf')

        # merge name
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_collectf')

        # merge uniprot_accession
        df = self.merge_columns(df=df, column='uniprot_accession',
                                left='uniprot_accession_annotation', right='uniprot_accession_collectf')

        # merge ncbi_protein
        df = self.merge_columns(df=df, column='ncbi_protein',
                                left='ncbi_protein_annotation', right='ncbi_protein_collectf')

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
