import pandas as pd

from protrend.model import Organism, Regulator, Gene
from protrend.report import ProtrendReporter
from protrend.transform.abasy.base import AbasyTransformer, AbasyConnector
from protrend.transform.mix_ins import OrganismMixIn
from protrend.utils import SetList
from protrend.utils.processors import to_int_str


class OrganismTransformer(OrganismMixIn, AbasyTransformer,
                          source='abasy',
                          version='0.0.0',
                          node=Organism,
                          order=100,
                          register=True):
    species = ('Bacillus subtilis',
               'Escherichia coli',
               'Corynebacterium glutamicum',
               'Mycobacterium tuberculosis',
               'Pseudomonas aeruginosa')

    strain = ('subsp. subtilis str. 168',
              'str. K-12 substr. MG1655',
              'ATCC 13032',
              'H37Rv',
              'PAO1')

    ncbi_taxonomy = (224308,
                     511145,
                     196627,
                     83332,
                     208964)

    refseq_accession = ('GCF_000009045.1',
                        'GCF_000005845.2',
                        'GCF_000196335.1',
                        'GCF_000195955.2',
                        'GCF_000006765.1')

    refseq_ftp = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/196/335/GCF_000196335.1_ASM19633v1',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1')

    genbank_accession = ('GCA_000009045.1',
                         'GCA_000005845.2',
                         'GCA_000196335.1',
                         'GCA_000195955.2',
                         'GCA_000006765.1')

    genbank_ftp = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/045/GCA_000009045.1_ASM904v1',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/196/335/GCA_000196335.1_ASM19633v1',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/765/GCA_000006765.1_ASM676v1')

    ncbi_assembly = (30588,
                     79781,
                     264778,
                     538048,
                     28348)

    assembly_accession = ('ASM904v1',
                          'ASM584v2',
                          'ASM19633v1',
                          'ASM19595v2',
                          'ASM676v1')

    name = ('Bacillus subtilis subsp. subtilis str. 168',
            'Escherichia coli str. K-12 substr. MG1655',
            'Corynebacterium glutamicum ATCC 13032',
            'Mycobacterium tuberculosis H37Rv',
            'Pseudomonas aeruginosa PAO1')

    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession'])

    def transform(self):
        # 'bsub': '224308',
        #  'cglu': '196627',
        #  'ecol': '511145',
        #  'mtub': '83332',
        #  'paer': '208964'
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

        # 'saur': '367830',
        # 'scoe': '100226',
        # 'spne': '406558',
        # 'spyo': '186103'
        organisms_to_annotate = dict(input_value=['367830', '100226', '406558', '186103'],
                                     ncbi_taxonomy=['367830', '100226', '406558', '186103'],
                                     name=['Staphylococcus aureus subsp. aureus USA300',
                                           'Streptomyces coelicolor A3(2)',
                                           'Streptococcus pneumoniae SP9-BS68',
                                           'Streptococcus pyogenes MGAS8232'])
        organisms_to_annotate = pd.DataFrame(organisms_to_annotate)

        organisms = self.annotate_organisms(organisms_to_annotate)
        organisms = organisms.drop(columns=['input_value'])

        df = pd.concat([annotated_organisms, organisms])

        ProtrendReporter.report_objects(source=self.source, version=self.version,
                                        system='extract', label=self.node.node_name(),
                                        objects=df.shape[0], properties=df.shape[1])

        self.stack_transformed_nodes(df)
        return df


class OrganismToRegulatorConnector(AbasyConnector,
                                   source='abasy',
                                   version='0.0.0',
                                   from_node=Organism,
                                   to_node=Regulator,
                                   register=True):

    def connect(self):
        source_processors = {'ncbi_taxonomy': [to_int_str]}
        target_processors = {'taxonomy': [to_int_str]}
        df = self.create_connection(source='organism', target='regulator',
                                    source_on='ncbi_taxonomy', target_on='taxonomy',
                                    source_processors=source_processors,
                                    target_processors=target_processors)
        self.stack_connections(df)


class OrganismToGeneConnector(AbasyConnector,
                              source='abasy',
                              version='0.0.0',
                              from_node=Organism,
                              to_node=Gene,
                              register=True):

    def connect(self):
        source_processors = {'ncbi_taxonomy': [to_int_str]}
        target_processors = {'taxonomy': [to_int_str]}
        df = self.create_connection(source='organism', target='gene',
                                    source_on='ncbi_taxonomy', target_on='taxonomy',
                                    source_processors=source_processors,
                                    target_processors=target_processors)
        self.stack_connections(df)
