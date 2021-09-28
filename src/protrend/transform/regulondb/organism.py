import pandas as pd

from protrend.model.model import Source
from protrend.transform.regulondb.base import RegulondbTransformer
from protrend.utils import SetList


class OrganismTransformer(RegulondbTransformer):
    species = 'Escherichia coli'
    strain = 'str. K-12 substr. MG1655'
    ncbi_taxonomy = 511145
    refseq_accession = 'GCF_000005845.2'
    refseq_ftp = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2'
    genbank_accession = 'GCA_000005845.2'
    genbank_ftp = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2'
    ncbi_assembly = 79781
    assembly_accession = 'ASM584v2'
    name = 'Escherichia coli str. K-12 substr. MG1655'

    default_node = Source
    default_order = 100
    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession'])

    def transform(self):
        ecoli = dict(name=[self.name],
                     species=[self.species],
                     strain=[self.strain],
                     ncbi_taxonomy=[self.ncbi_taxonomy],
                     refseq_accession=[self.refseq_accession],
                     refseq_ftp=[self.refseq_ftp],
                     genbank_accession=[self.genbank_accession],
                     genbank_ftp=[self.genbank_ftp],
                     ncbi_assembly=[self.ncbi_assembly],
                     assembly_accession=[self.assembly_accession])

        df = pd.DataFrame(ecoli, index=[0])

        self._stack_transformed_nodes(df)

        return df
