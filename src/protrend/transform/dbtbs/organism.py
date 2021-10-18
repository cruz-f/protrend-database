import pandas as pd

from protrend.model.model import Organism
from protrend.transform.dbtbs.base import DBTBSTransformer
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
