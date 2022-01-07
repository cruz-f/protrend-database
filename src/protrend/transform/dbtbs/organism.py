from protrend.transform import BaseOrganismTransformer
from protrend.model import Organism, Regulator, Gene, TFBS, RegulatoryInteraction
from protrend.transform.dbtbs.base import DBTBSTransformer, DBTBSConnector
from protrend.utils import SetList


class OrganismTransformer(DBTBSTransformer, BaseOrganismTransformer,
                          source='dbtbs',
                          version='0.0.4',
                          node=Organism,
                          order=100,
                          register=True):
    species = ['Bacillus subtilis']
    strain = ['subsp. subtilis str. 168']
    ncbi_taxonomy = [224308]
    refseq_accession = ['GCF_000009045.1']
    refseq_ftp = ['ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1']
    genbank_accession = ['GCA_000009045.1']
    genbank_ftp = ['ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/045/GCA_000009045.1_ASM904v1']
    ncbi_assembly = [30588]
    assembly_accession = ['ASM904v1']
    name = ['Bacillus subtilis subsp. subtilis str. 168']

    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession'])


class OrganismToRegulatorConnector(DBTBSConnector,
                                   source='dbtbs',
                                   version='0.0.4',
                                   from_node=Organism,
                                   to_node=Regulator,
                                   register=True):
    default_connect_stack = {'organism': 'integrated_organism.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        df = self.create_connection(source='organism', target='regulator', cardinality='one_to_many')
        self.stack_json(df)


class OrganismToGeneConnector(DBTBSConnector,
                              source='dbtbs',
                              version='0.0.4',
                              from_node=Organism,
                              to_node=Gene,
                              register=True):
    default_connect_stack = {'organism': 'integrated_organism.json', 'gene': 'integrated_gene.json'}

    def connect(self):
        df = self.create_connection(source='organism', target='gene', cardinality='one_to_many')
        self.stack_json(df)


class OrganismToTFBSConnector(DBTBSConnector,
                              source='dbtbs',
                              version='0.0.4',
                              from_node=Organism,
                              to_node=TFBS,
                              register=True):
    default_connect_stack = {'organism': 'integrated_organism.json', 'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        df = self.create_connection(source='organism', target='tfbs', cardinality='one_to_many')
        self.stack_json(df)


class OrganismToRegulatoryInteractionConnector(DBTBSConnector,
                                               source='dbtbs',
                                               version='0.0.4',
                                               from_node=Organism,
                                               to_node=RegulatoryInteraction,
                                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin', source_column='organism')
        self.stack_json(df)
