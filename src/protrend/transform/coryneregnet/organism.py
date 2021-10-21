import pandas as pd

from protrend.model.model import Organism
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer
from protrend.utils import SetList


class OrganismTransformer(CoryneRegNetTransformer):
    species = ('Bacillus subtilis',
               'Escherichia coli',
               'Corynebacterium glutamicum',
               'Mycobacterium tuberculosis')

    strain = ('subsp. subtilis str. 168',
              'str. K-12 substr. MG1655',
              'ATCC 13032',
              'H37Rv')

    ncbi_taxonomy = (224308,
                     511145,
                     196627,
                     83332)

    refseq_accession = ('GCF_000009045.1',
                        'GCF_000005845.2',
                        'GCF_000196335.1',
                        'GCF_000195955.2')

    refseq_ftp = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/196/335/GCF_000196335.1_ASM19633v1',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2')

    genbank_accession = ('GCA_000009045.1',
                         'GCA_000005845.2',
                         'GCA_000196335.1',
                         'GCA_000195955.2')

    genbank_ftp = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/045/GCA_000009045.1_ASM904v1',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/196/335/GCA_000196335.1_ASM19633v1',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2')

    ncbi_assembly = (30588,
                     79781,
                     264778,
                     538048)

    assembly_accession = ('ASM904v1',
                          'ASM584v2',
                          'ASM19633v1',
                          'ASM19595v2')

    name = ('Bacillus subtilis subsp. subtilis str. 168',
            'Escherichia coli str. K-12 substr. MG1655',
            'Corynebacterium glutamicum ATCC 13032',
            'Mycobacterium tuberculosis H37Rv')

    default_node = Organism
    default_order = 100
    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession'])

    def transform(self):
        org = dict(name=self.name,
                   species=self.species,
                   strain=self.strain,
                   ncbi_taxonomy=self.ncbi_taxonomy,
                   refseq_accession=self.refseq_accession,
                   refseq_ftp=self.refseq_ftp,
                   genbank_accession=self.genbank_accession,
                   genbank_ftp=self.genbank_ftp,
                   ncbi_assembly=self.ncbi_assembly,
                   assembly_accession=self.assembly_accession)

        df = pd.DataFrame(org)

        self._stack_transformed_nodes(df)

        return df