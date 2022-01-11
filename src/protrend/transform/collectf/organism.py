from typing import List

import pandas as pd

from protrend.bioapis import entrez_summary
from protrend.io import read_json_lines, read_from_stack
from protrend.model import Organism, Regulator, Gene, TFBS, RegulatoryInteraction
from protrend.transform.collectf.base import CollectfTransformer, CollectfConnector
from protrend.transform.mix_ins import OrganismMixIn
from protrend.transform.transformations import (merge_columns, create_input_value, drop_duplicates, group_by,
                                                drop_empty_string)
from protrend.utils import SetList, is_null
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_int_str, take_last, flatten_set_list


class OrganismTransformer(OrganismMixIn, CollectfTransformer,
                          source='collectf',
                          version='0.0.1',
                          node=Organism,
                          order=100,
                          register=True):
    default_transform_stack = {'organism': 'Organism.json'}
    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession',
                       'name_collectf', 'genome_accession', 'taxonomy', 'regulon', 'tfbs'])
    read_columns = SetList(['name', 'genome_accession', 'taxonomy', 'regulon', 'tfbs'])

    @staticmethod
    def get_ncbi_taxa_from_nucleotide(nucleotide: List[str]):
        ncbi_taxa = []
        for acc in nucleotide:

            ncbi_taxon = None

            if not is_null(acc):
                summary = entrez_summary(identifier=acc, db='nucleotide')

                if 'TaxId' in summary:
                    ncbi_taxon = str(summary['TaxId'])

            ncbi_taxa.append(ncbi_taxon)

        return ncbi_taxa

    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        organism = apply_processors(organism, name=[rstrip, lstrip])
        organism = organism.dropna(subset=['name'])
        organism = drop_empty_string(organism, 'name')

        aggregation = {'genome_accession': take_last, 'taxonomy': take_last}
        organism = group_by(df=organism, column='name', aggregation=aggregation, default=flatten_set_list)

        organism = apply_processors(organism, genome_accession=[rstrip, lstrip], taxonomy=to_int_str)
        organism = drop_duplicates(df=organism, subset=['genome_accession', 'name'])

        nucleotide = organism['genome_accession'].to_list()
        ncbi_taxa = self.get_ncbi_taxa_from_nucleotide(nucleotide)
        organism = organism.assign(ncbi_taxonomy=ncbi_taxa)

        organism = create_input_value(organism, col='name')
        return organism

    def transform(self):
        organism = read_from_stack(stack=self.transform_stack, key='organism',
                                   columns=self.read_columns, reader=read_json_lines)

        organisms = self.transform_organism(organism)
        annotated_organisms = self.annotate_organisms(organisms)

        df = pd.merge(annotated_organisms, organisms, on='input_value', suffixes=('_annotation', '_collectf'))

        # merge name
        df = df.assign(name_to_merge=df['name_collectf'])
        df = merge_columns(df=df, column='name', left='name_annotation', right='name_to_merge')

        df = apply_processors(df, ncbi_taxonomy=to_int_str, ncbi_assembly=to_int_str)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df


class OrganismToRegulatorConnector(CollectfConnector,
                                   source='collectf',
                                   version='0.0.1',
                                   from_node=Organism,
                                   to_node=Regulator,
                                   register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='organism', target_column='regulator')
        self.stack_json(df)


class OrganismToGeneConnector(CollectfConnector,
                              source='collectf',
                              version='0.0.1',
                              from_node=Organism,
                              to_node=Gene,
                              register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='organism', target_column='gene')
        self.stack_json(df)


class OrganismToTFBSConnector(CollectfConnector,
                              source='collectf',
                              version='0.0.1',
                              from_node=Organism,
                              to_node=TFBS,
                              register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='organism', target_column='tfbs')
        self.stack_json(df)


class OrganismToRegulatoryInteractionConnector(CollectfConnector,
                                               source='collectf',
                                               version='0.0.1',
                                               from_node=Organism,
                                               to_node=RegulatoryInteraction,
                                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='organism')
        self.stack_json(df)
