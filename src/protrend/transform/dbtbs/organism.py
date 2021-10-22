import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Organism, Regulator, Operon, Gene, TFBS, RegulatoryInteraction
from protrend.transform.dbtbs.base import DBTBSTransformer, DBTBSConnector
from protrend.transform.dbtbs.gene import GeneTransformer
from protrend.transform.dbtbs.operon import OperonTransformer
from protrend.transform.dbtbs.regulator import RegulatorTransformer
from protrend.transform.dbtbs.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.transform.dbtbs.tfbs import TFBSTransformer
from protrend.utils import SetList


class OrganismTransformer(DBTBSTransformer):
    species = 'Bacillus subtilis'
    strain = 'subsp. subtilis str. 168'
    ncbi_taxonomy = 224308
    refseq_accession = 'GCF_000009045.1'
    refseq_ftp = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1'
    genbank_accession = 'GCA_000009045.1'
    genbank_ftp = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/045/GCA_000009045.1_ASM904v1'
    ncbi_assembly = 30588
    assembly_accession = 'ASM904v1'
    name = 'Bacillus subtilis subsp. subtilis str. 168'

    default_node = Organism
    default_order = 100
    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession'])

    def transform(self):
        org = dict(name=[self.name],
                   species=[self.species],
                   strain=[self.strain],
                   ncbi_taxonomy=[self.ncbi_taxonomy],
                   refseq_accession=[self.refseq_accession],
                   refseq_ftp=[self.refseq_ftp],
                   genbank_accession=[self.genbank_accession],
                   genbank_ftp=[self.genbank_ftp],
                   ncbi_assembly=[self.ncbi_assembly],
                   assembly_accession=[self.assembly_accession])

        df = pd.DataFrame(org, index=[0])

        self._stack_transformed_nodes(df)

        return df


class RegulatorToOrganismConnector(DBTBSConnector):
    default_from_node = Regulator
    default_to_node = Organism
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)

        from_identifiers = regulator['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = organism['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class OperonToOrganismConnector(DBTBSConnector):
    default_from_node = Operon
    default_to_node = Organism
    default_connect_stack = {'operon': 'integrated_operon.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)

        from_identifiers = operon['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = organism['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class GeneToOrganismConnector(DBTBSConnector):
    default_from_node = Gene
    default_to_node = Organism
    default_connect_stack = {'gene': 'integrated_gene.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)

        from_identifiers = gene['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = organism['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class TFBSToOrganismConnector(DBTBSConnector):
    default_from_node = TFBS
    default_to_node = Organism
    default_connect_stack = {'tfbs': 'integrated_tfbs.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)

        from_identifiers = tfbs['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = organism['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToOrganismConnector(DBTBSConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Organism
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json',
                             'organism': 'integrated_organism.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = organism['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
