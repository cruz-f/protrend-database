from typing import List, Union

import pandas as pd
from Bio.SeqRecord import SeqRecord

from protrend.bioapis import map_uniprot_identifiers, fetch_uniprot_record
from protrend.io import read_from_stack, read_json_lines, read_json_frame
from protrend.model.model import Regulator
from protrend.transform import GeneDTO
from protrend.transform.annotation import annotate_genes
from protrend.transform.collectf.base import CollectfTransformer
from protrend.transform.processors import take_first, flatten_set, apply_processors, to_list_nan


def _find_in_mapping(acc: str, mapping: pd.DataFrame) -> List[Union[str, None]]:
    acc_mask = mapping['From'] == acc
    acc_series: pd.Series = mapping.loc[acc_mask, 'To']

    res = list(acc_series)

    if res:
        return res

    return [None]


class RegulatorTransformer(CollectfTransformer):
    default_node = Regulator
    default_node_factors = ('uniprot_accession', 'ncbi_protein', 'ncbi_gene',
                            'genbank_accession', 'refseq_accession',
                            'locus_tag')
    default_transform_stack = {'regulon': 'Regulon.json', 'organism': 'integrated_organism.json'}
    default_order = 90
    columns = {'protrend_id',
               'mechanism', 'name',
               'locus_tag', 'synonyms', 'function', 'description',
               'ncbi_gene', 'ncbi_protein',
               'genbank_accession', 'refseq_accession',
               'sequence', 'strand', 'start', 'stop',
               'uniprot_accession', 'name', 'url', 'organism', 'operon', 'gene', 'tfbs', 'experimental_evidence',
               'organism_protrend_id', 'organism_name_collectf', 'ncbi_taxonomy'}
    read_columns = {'uniprot_accession', 'name', 'url', 'organism', 'operon', 'gene', 'tfbs', 'experimental_evidence'}

    def _transform_regulon(self, regulon: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
        regulon = apply_processors(regulon, tfbs=to_list_nan, experimental_evidence=to_list_nan,
                                   operon=to_list_nan, gene=to_list_nan)
        aggregation = {'tfbs': flatten_set, 'experimental_evidence': flatten_set,
                       'operon': flatten_set, 'gene': flatten_set}
        regulon = self.group_by(df=regulon, column='uniprot_accession', aggregation=aggregation, default=take_first)

        regulon['mechanism'] = 'transcription factor'

        df = pd.merge(regulon, organism, how='left', left_on='organism', right_on='organism_name_collectf')

        df = self.drop_duplicates(df=df, subset=['uniprot_accession', 'organism'],
                                  perfect_match=True, preserve_nan=True)
        df = self.create_input_value(df=df, col='uniprot_accession')
        return df

    @staticmethod
    def _locus_tag_from_uniprot_record(record: SeqRecord) -> Union[str, None]:
        annotations = getattr(record, 'annotations', {})

        if annotations:
            loci = annotations.get('gene_name_ordered locus', [])

            if not loci:
                loci = annotations.get('gene_name_ORF', [])

            for locus in loci:
                return locus

        return None

    def _annotate_tfs(self,
                      accessions: List[Union[None, str]],
                      names: List[str],
                      taxa: List[str]) -> pd.DataFrame:

        # map uniprot_accessions to ncbi_proteins
        uniprot_ncbi_proteins = map_uniprot_identifiers(accessions, from_='ACC', to='P_GI')
        ncbi_proteins = [_find_in_mapping(accession, uniprot_ncbi_proteins)[0] for accession in accessions]

        # fetch uniprot record to retrieve locus_tag
        uniprot_records = [fetch_uniprot_record(accession) for accession in accessions]
        loci = [self._locus_tag_from_uniprot_record(record) for record in uniprot_records]

        # annotate with uniprot_accessions, ncbi_proteins, locus_tag for ncbi_gene
        dtos = [GeneDTO(input_value=accession) for accession in accessions]
        annotate_genes(dtos=dtos, loci=loci, names=names, taxa=taxa,
                       uniprot_proteins=accessions, ncbi_proteins=ncbi_proteins)

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

    def transform(self):
        regulon = read_from_stack(stack=self.transform_stack, file='regulon',
                                  default_columns=self.read_columns, reader=read_json_lines)

        from protrend.transform.collectf import OrganismTransformer
        organism = read_from_stack(stack=self.transform_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = self.select_columns(organism, 'protrend_id', 'name_collectf', 'ncbi_taxonomy')
        organism = organism.rename(columns={'protrend_id': 'organism_protrend_id',
                                            'name_collectf': 'organism_name_collectf'})

        regulon = self._transform_regulon(regulon, organism)

        accessions = regulon['uniprot_accession'].tolist()
        names = regulon['name'].tolist()
        taxa = regulon['ncbi_taxonomy'].tolist()
        tfs = self._annotate_tfs(accessions=accessions, names=names, taxa=taxa)

        tf = pd.merge(tfs, regulon, on='input_value', suffixes=('_annotation', '_collectf'))

        tf = self.merge_columns(df=tf, column='name', left='name_annotation', right='name_collectf')
        tf = self.merge_columns(df=tf, column='uniprot_accession',
                                left='uniprot_accession_annotation', right='uniprot_accession_collectf')

        # dropping those tfs with incomplete locus_tag, as it means that the uniprot accession is wrongly assigned by
        # collectf or absent
        tf = tf.dropna(subset=['locus_tag'])
        non_empty_tf = tf['locus_tag'] != ''
        tf = tf[non_empty_tf]

        tf = tf.drop(columns=['input_value'])

        self._stack_transformed_nodes(tf)
        return tf
