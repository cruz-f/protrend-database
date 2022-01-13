from protrend.model import Regulator, Organism, Gene, TFBS, RegulatoryInteraction
from protrend.transform.mix_ins import OrganismMixIn
from protrend.transform.regulondb.base import RegulonDBTransformer, RegulonDBConnector
from protrend.utils import SetList


class OrganismTransformer(OrganismMixIn, RegulonDBTransformer,
                          source='regulondb',
                          version='0.0.0',
                          node=Organism,
                          order=100,
                          register=True):
    species = ['Escherichia coli']
    strain = ['str. K-12 substr. MG1655']
    ncbi_taxonomy = [511145]
    refseq_accession = ['GCF_000005845.2']
    refseq_ftp = ['ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2']
    genbank_accession = ['GCA_000005845.2']
    genbank_ftp = ['ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2']
    ncbi_assembly = [79781]
    assembly_accession = ['ASM584v2']
    name = ['Escherichia coli str. K-12 substr. MG1655']

    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession'])


class OrganismToRegulatorConnector(RegulonDBConnector,
                                   source='regulondb',
                                   version='0.0.0',
                                   from_node=Organism,
                                   to_node=Regulator,
                                   register=True):
    default_connect_stack = {'organism': 'integrated_organism.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        df = self.create_connection(source='organism', target='regulator', cardinality='one_to_many')
        self.stack_connections(df)


class OrganismToGeneConnector(RegulonDBConnector,
                              source='regulondb',
                              version='0.0.0',
                              from_node=Organism,
                              to_node=Gene,
                              register=True):
    default_connect_stack = {'organism': 'integrated_organism.json', 'gene': 'integrated_gene.json'}

    def connect(self):
        df = self.create_connection(source='organism', target='gene', cardinality='one_to_many')
        self.stack_connections(df)


class OrganismToTFBSConnector(RegulonDBConnector,
                              source='regulondb',
                              version='0.0.0',
                              from_node=Organism,
                              to_node=TFBS,
                              register=True):
    default_connect_stack = {'organism': 'integrated_organism.json', 'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        df = self.create_connection(source='organism', target='tfbs', cardinality='one_to_many')
        self.stack_connections(df)


class OrganismToRegulatoryInteractionConnector(RegulonDBConnector,
                                               source='regulondb',
                                               version='0.0.0',
                                               from_node=Organism,
                                               to_node=RegulatoryInteraction,
                                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin', source_column='organism')
        self.stack_connections(df)
