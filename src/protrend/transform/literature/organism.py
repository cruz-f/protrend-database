from protrend.model import Organism, Regulator, Gene, RegulatoryInteraction
from protrend.transform.literature.base import LiteratureTransformer, LiteratureConnector
from protrend.transform.mix_ins import OrganismMixIn
from protrend.utils import SetList


class OrganismTransformer(OrganismMixIn, LiteratureTransformer,
                          source='literature',
                          version='0.0.0',
                          node=Organism,
                          order=100,
                          register=True):
    species = ['Bacillus subtilis',
               'Escherichia coli',
               'Mycobacterium tuberculosis']

    strain = ['subsp. subtilis str. 168',
              'str. K-12 substr. MG1655',
              'H37Rv']

    ncbi_taxonomy = [224308,
                     511145,
                     83332]

    refseq_accession = ['GCF_000009045.1',
                        'GCF_000005845.2',
                        'GCF_000195955.2']

    refseq_ftp = ['ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2']

    genbank_accession = ['GCA_000009045.1',
                         'GCA_000005845.2',
                         'GCA_000195955.2']

    genbank_ftp = ['ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/045/GCA_000009045.1_ASM904v1',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2']

    ncbi_assembly = [30588,
                     79781,
                     538048]

    assembly_accession = ['ASM904v1',
                          'ASM584v2',
                          'ASM19595v2']

    name = ['Bacillus subtilis subsp. subtilis str. 168',
            'Escherichia coli str. K-12 substr. MG1655',
            'Mycobacterium tuberculosis H37Rv']

    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession'])


class OrganismToRegulatorConnector(LiteratureConnector,
                                   source='literature',
                                   version='0.0.0',
                                   from_node=Organism,
                                   to_node=Regulator,
                                   register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='organism', target_column='regulator')
        self.stack_connections(df)


class OrganismToGeneConnector(LiteratureConnector,
                              source='literature',
                              version='0.0.0',
                              from_node=Organism,
                              to_node=Gene,
                              register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='organism', target_column='gene')
        self.stack_connections(df)


class OrganismToRegulatoryInteractionConnector(LiteratureConnector,
                                               source='literature',
                                               version='0.0.0',
                                               from_node=Organism,
                                               to_node=RegulatoryInteraction,
                                               register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='organism')
        self.stack_connections(df)
