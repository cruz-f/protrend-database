from typing import List

import pandas as pd

from protrend.bioapis import entrez_summary
from protrend.io import read_json_lines, read_from_stack, read_json_frame
from protrend.model.model import Organism, Regulator, Operon, Gene, TFBS, RegulatoryInteraction
from protrend.transform import OrganismDTO
from protrend.transform.annotation import annotate_organisms
from protrend.transform.collectf.base import CollectfTransformer, CollectfConnector
from protrend.transform.collectf.gene import GeneTransformer
from protrend.transform.collectf.operon import OperonTransformer
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.transform.collectf.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.transform.collectf.tfbs import TFBSTransformer
from protrend.transform.processors import apply_processors, rstrip, lstrip, to_int_str, take_last, flatten_set, to_list
from protrend.utils.miscellaneous import is_null


class OrganismTransformer(CollectfTransformer):
    default_node = Organism
    default_node_factors = ('ncbi_taxonomy', 'name')
    default_transform_stack = {'organism': 'Organism.json'}
    default_order = 100
    columns = {'protrend_id',
               'genome_accession', 'taxonomy', 'regulon', 'tfbs', 'name_collectf'
                                                                  'name', 'species', 'strain',
               'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
               'genbank_accession', 'genbank_ftp',
               'ncbi_assembly', 'assembly_accession'}
    read_columns = {'name', 'genome_accession', 'taxonomy', 'regulon', 'tfbs'}

    def _transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        aggregation = {'genome_accession': take_last, 'taxonomy': take_last}
        organism = self.group_by(df=organism, column='name', aggregation=aggregation, default=flatten_set)

        organism = apply_processors(organism, name=[rstrip, lstrip], genome_accession=[rstrip, lstrip],
                                    taxonomy=to_int_str)

        organism = self.create_input_value(organism, col='name')

        return organism

    @staticmethod
    def _transform_organisms(nucleotide: List[str], names: List[str]):

        identifiers = []
        for acc in nucleotide:

            ncbi_tax = None

            if not is_null(acc):
                summary = entrez_summary(identifier=acc, db='nucleotide')

                if 'TaxId' in summary:
                    ncbi_tax = str(summary['TaxId'])

            identifiers.append(ncbi_tax)

        dtos = [OrganismDTO(input_value=name) for name in names]
        annotate_organisms(dtos=dtos, identifiers=identifiers, names=names)

        # name: List[str]
        # species: List[str]
        # strain: List[str]
        # ncbi_taxonomy: List[int]
        # refseq_accession: List[str]
        # refseq_ftp: List[str]
        # genbank_accession: List[str]
        # genbank_ftp: List[str]
        # ncbi_assembly: List[str]
        # assembly_accession: List[str]

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        organism = read_from_stack(stack=self.transform_stack, file='organism',
                                   default_columns=self.read_columns, reader=read_json_lines)
        organism = self._transform_organism(organism)

        names = organism['input_value'].tolist()
        nucleotide = organism['genome_accession'].tolist()
        organisms = self._transform_organisms(nucleotide, names)

        df = pd.merge(organisms, organism, on='input_value', suffixes=('_annotation', '_collectf'))

        names_collectf = df['name_collectf'].tolist()
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_collectf')
        df['name_collectf'] = names_collectf

        df = df.drop(columns=['input_value'])

        df = apply_processors(df, ncbi_taxonomy=to_int_str, ncbi_assembly=to_int_str)

        self._stack_transformed_nodes(df)

        return df


class OrganismToRegulatorConnector(CollectfConnector):
    default_from_node = Organism
    default_to_node = Regulator
    default_connect_stack = {'regulator': 'integrated_regulator.json'}

    def connect(self):
        regulator = read_from_stack(stack=self.connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        regulator = regulator.dropna(subset=['organism_protrend_id'])
        regulator = regulator.dropna(subset=['regulator_protrend_id'])
        regulator = regulator.drop_duplicates(subset=['organism_protrend_id', 'regulator_protrend_id'])

        from_identifiers = regulator['organism_protrend_id'].tolist()
        to_identifiers = regulator['regulator_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)


class OrganismToOperonConnector(CollectfConnector):
    default_from_node = Organism
    default_to_node = Operon
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'operon': 'integrated_operon.json'}

    def connect(self):
        regulator = read_from_stack(stack=self.connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        regulator = regulator.dropna(subset=['organism_protrend_id'])
        regulator = regulator.dropna(subset=['regulator_protrend_id'])
        regulator = regulator.drop_duplicates(subset=['organism_protrend_id', 'regulator_protrend_id'])

        operon = read_from_stack(stack=self.connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = operon.rename(columns={'protrend_id': 'operon_protrend_id'})
        operon = apply_processors(operon, regulon=to_list)
        operon = operon.explode(column='regulon')

        df = pd.merge(regulator, operon, left_on='uniprot_accession', right_on='regulon')
        df = df.dropna(subset=['organism_protrend_id'])
        df = df.dropna(subset=['operon_protrend_id'])
        df = df.drop_duplicates(subset=['organism_protrend_id', 'operon_protrend_id'])

        from_identifiers = df['organism_protrend_id'].tolist()
        to_identifiers = df['operon_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)


class OrganismToGeneConnector(CollectfConnector):
    default_from_node = Organism
    default_to_node = Gene
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'gene': 'integrated_gene.json'}

    def connect(self):
        regulator = read_from_stack(stack=self.connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        regulator = regulator.dropna(subset=['organism_protrend_id'])
        regulator = regulator.dropna(subset=['regulator_protrend_id'])
        regulator = regulator.drop_duplicates(subset=['organism_protrend_id', 'regulator_protrend_id'])

        gene = read_from_stack(stack=self.connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id'})
        gene = apply_processors(gene, regulon=to_list)
        gene = gene.explode(column='regulon')

        df = pd.merge(regulator, gene, left_on='uniprot_accession', right_on='regulon')
        df = df.dropna(subset=['organism_protrend_id'])
        df = df.dropna(subset=['gene_protrend_id'])
        df = df.drop_duplicates(subset=['organism_protrend_id', 'gene_protrend_id'])

        from_identifiers = df['organism_protrend_id'].tolist()
        to_identifiers = df['gene_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)


class OrganismToTFBSConnector(CollectfConnector):
    default_from_node = Organism
    default_to_node = TFBS
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        regulator = read_from_stack(stack=self.connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        regulator = regulator.dropna(subset=['organism_protrend_id'])
        regulator = regulator.dropna(subset=['regulator_protrend_id'])
        regulator = regulator.drop_duplicates(subset=['organism_protrend_id', 'regulator_protrend_id'])

        tfbs = read_from_stack(stack=self.connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs_protrend_id'})
        tfbs = apply_processors(tfbs, regulon=to_list)
        tfbs = tfbs.explode(column='regulon')

        df = pd.merge(regulator, tfbs, left_on='uniprot_accession', right_on='regulon')
        df = df.dropna(subset=['organism_protrend_id'])
        df = df.dropna(subset=['tfbs_protrend_id'])
        df = df.drop_duplicates(subset=['organism_protrend_id', 'tfbs_protrend_id'])

        from_identifiers = df['organism_protrend_id'].tolist()
        to_identifiers = df['tfbs_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)


class OrganismToRegulatoryInteractionConnector(CollectfConnector):
    default_from_node = Organism
    default_to_node = RegulatoryInteraction
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'rin': 'integrated_tfbs.json'}

    def connect(self):
        regulator = read_from_stack(stack=self.connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        regulator = regulator.dropna(subset=['organism_protrend_id'])
        regulator = regulator.dropna(subset=['regulator_protrend_id'])
        regulator = regulator.drop_duplicates(subset=['organism_protrend_id', 'regulator_protrend_id'])

        rin = read_from_stack(stack=self.connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        rin = rin.rename(columns={'protrend_id': 'rin_protrend_id'})
        rin = apply_processors(rin, regulon=to_list)
        rin = rin.explode(column='regulon')

        df = pd.merge(regulator, rin, left_on='uniprot_accession', right_on='regulon')
        df = df.dropna(subset=['organism_protrend_id'])
        df = df.dropna(subset=['rin_protrend_id'])
        df = df.drop_duplicates(subset=['organism_protrend_id', 'rin_protrend_id'])

        from_identifiers = df['organism_protrend_id'].tolist()
        to_identifiers = df['rin_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)
