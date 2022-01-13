import pandas as pd

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
               'Mycobacterium tuberculosis',
               'Pseudomonas aeruginosa']

    strain = ['subsp. subtilis str. 168',
              'str. K-12 substr. MG1655',
              'H37Rv',
              'PAO1']

    ncbi_taxonomy = [224308,
                     511145,
                     83332,
                     208964]

    refseq_accession = ['GCF_000009045.1',
                        'GCF_000005845.2',
                        'GCF_000195955.2',
                        'GCF_000006765.1']

    refseq_ftp = ['ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1']

    genbank_accession = ['GCA_000009045.1',
                         'GCA_000005845.2',
                         'GCA_000195955.2',
                         'GCA_000006765.1']

    genbank_ftp = ['ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/045/GCA_000009045.1_ASM904v1',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/765/GCA_000006765.1_ASM676v1']

    ncbi_assembly = [30588,
                     79781,
                     538048,
                     28348]

    assembly_accession = ['ASM904v1',
                          'ASM584v2',
                          'ASM19595v2',
                          'ASM676v1']

    name = ['Bacillus subtilis subsp. subtilis str. 168',
            'Escherichia coli str. K-12 substr. MG1655',
            'Mycobacterium tuberculosis H37Rv',
            'Pseudomonas aeruginosa PAO1']

    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession'])

    def transform(self):
        # 'bsub': '224308',
        #  'ecol': '511145',
        #  'mtub': '83332',
        #  'paer_PAO1': '208964'
        annotated_organisms = dict(name=self.name,
                                   species=self.species,
                                   strain=self.strain,
                                   ncbi_taxonomy=self.ncbi_taxonomy,
                                   refseq_accession=self.refseq_accession,
                                   refseq_ftp=self.refseq_ftp,
                                   genbank_accession=self.genbank_accession,
                                   genbank_ftp=self.genbank_ftp,
                                   ncbi_assembly=self.ncbi_assembly,
                                   assembly_accession=self.assembly_accession)

        annotated_organisms = pd.DataFrame(annotated_organisms)

        # 'paer_PAK': '1009714',
        # 'paer_PA14': '652611',
        # 'paer_PA103': '1081927',
        organisms_to_annotate = dict(input_values=['1009714', '652611', '1081927'],
                                     ncbi_taxonomy=['1009714', '652611', '1081927'],
                                     name=['Pseudomonas aeruginosa PAK', 'Pseudomonas aeruginosa PA14',
                                           'Pseudomonas aeruginosa PA103'])
        organisms_to_annotate = pd.DataFrame(organisms_to_annotate)

        organisms = self.annotate_organisms(organisms_to_annotate)
        organisms = organisms.drop(columns=['input_value'])

        df = pd.concat([annotated_organisms, organisms])

        self.stack_transformed_nodes(df)
        return df


class OrganismToRegulatorConnector(LiteratureConnector,
                                   source='literature',
                                   version='0.0.0',
                                   from_node=Organism,
                                   to_node=Regulator,
                                   register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

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
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

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
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='organism')
        self.stack_connections(df)
