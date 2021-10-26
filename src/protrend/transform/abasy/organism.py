from typing import List

import pandas as pd

from protrend.model.model import Organism
from protrend.transform import OrganismDTO
from protrend.transform.abasy.base import AbasyTransformer
from protrend.transform.annotation import annotate_organisms
from protrend.utils import SetList


class OrganismTransformer(AbasyTransformer):
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

    default_node = Organism
    default_order = 100
    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession'])

    @staticmethod
    def _transform_organisms(identifiers: List[str], names: List[str]):
        dtos = [OrganismDTO(input_value=identifier) for identifier in identifiers]
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

        annotated_organisms_df = pd.DataFrame(annotated_organisms)

        # 'saur': '367830',
        # 'scoe': '100226',
        # 'spne': '406558',
        # 'spyo': '186103'
        tax_ids = ['367830', '100226', '406558', '186103']
        names = ['Staphylococcus aureus subsp. aureus USA300', 'Streptomyces coelicolor A3(2)',
                 'Streptococcus pneumoniae SP9-BS68', 'Streptococcus pyogenes MGAS8232']
        organisms = self._transform_organisms(identifiers=tax_ids, names=names)
        organisms = organisms.drop(columns=['input_value'])

        df = pd.concat([annotated_organisms_df, organisms], axis=0)

        self._stack_transformed_nodes(df)

        return df
