from typing import List

import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Organism, Regulator, Operon, Gene, RegulatoryInteraction, Effector
from protrend.transform import OrganismDTO
from protrend.annotation import annotate_organisms
from protrend.transform.literature.base import LiteratureTransformer, LiteratureConnector
from protrend.transform.literature.effector import EffectorTransformer
from protrend.transform.literature.operon import OperonTransformer
from protrend.transform.literature.regulator import RegulatorTransformer
from protrend.transform.literature.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.utils.processors import apply_processors, to_int_str
from protrend.utils import SetList


class OrganismTransformer(LiteratureTransformer):
    species = ('Bacillus subtilis',
               'Escherichia coli',
               'Mycobacterium tuberculosis',
               'Pseudomonas aeruginosa')

    strain = ('subsp. subtilis str. 168',
              'str. K-12 substr. MG1655',
              'H37Rv',
              'PAO1')

    ncbi_taxonomy = (224308,
                     511145,
                     83332,
                     208964)

    refseq_accession = ('GCF_000009045.1',
                        'GCF_000005845.2',
                        'GCF_000195955.2',
                        'GCF_000006765.1')

    refseq_ftp = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2',
                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1')

    genbank_accession = ('GCA_000009045.1',
                         'GCA_000005845.2',
                         'GCA_000195955.2',
                         'GCA_000006765.1')

    genbank_ftp = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/045/GCA_000009045.1_ASM904v1',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2',
                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/765/GCA_000006765.1_ASM676v1')

    ncbi_assembly = (30588,
                     79781,
                     538048,
                     28348)

    assembly_accession = ('ASM904v1',
                          'ASM584v2',
                          'ASM19595v2',
                          'ASM676v1')

    name = ('Bacillus subtilis subsp. subtilis str. 168',
            'Escherichia coli str. K-12 substr. MG1655',
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

        annotated_organisms_df = pd.DataFrame(annotated_organisms)

        # 'paer_PAK': '1009714',
        # 'paer_PA14': '652611',
        # 'paer_PA103': '1081927',
        tax_ids = ['1009714', '652611', '1081927']
        names = ['Pseudomonas aeruginosa PAK', 'Pseudomonas aeruginosa PA14', 'Pseudomonas aeruginosa PA103']
        organisms = self._transform_organisms(identifiers=tax_ids, names=names)
        organisms = organisms.drop(columns=['input_value'])

        df = pd.concat([annotated_organisms_df, organisms], axis=0)

        self._stack_transformed_nodes(df)

        return df


class RegulatorToOrganismConnector(LiteratureConnector):
    default_from_node = Regulator
    default_to_node = Organism
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = apply_processors(regulator, taxonomy=to_int_str)

        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = apply_processors(organism, ncbi_taxonomy=to_int_str)

        df = pd.merge(regulator, organism, left_on='taxonomy', right_on='ncbi_taxonomy',
                      suffixes=('_regulator', '_organism'))

        from_identifiers = df['protrend_id_regulator'].tolist()
        to_identifiers = df['protrend_id_organism'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class OperonToOrganismConnector(LiteratureConnector):
    default_from_node = Operon
    default_to_node = Organism
    default_connect_stack = {'operon': 'integrated_operon.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = apply_processors(operon, taxonomy=to_int_str)

        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = apply_processors(organism, ncbi_taxonomy=to_int_str)

        df = pd.merge(operon, organism, left_on='taxonomy', right_on='ncbi_taxonomy',
                      suffixes=('_operon', '_organism'))

        from_identifiers = df['protrend_id_operon'].tolist()
        to_identifiers = df['protrend_id_organism'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class GeneToOrganismConnector(LiteratureConnector):
    default_from_node = Gene
    default_to_node = Organism
    default_connect_stack = {'operon': 'integrated_operon.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = apply_processors(operon, taxonomy=to_int_str)

        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = apply_processors(organism, ncbi_taxonomy=to_int_str)

        df = pd.merge(operon, organism, left_on='taxonomy', right_on='ncbi_taxonomy',
                      suffixes=('_operon', '_organism'))
        df = df.explode(column='genes')

        from_identifiers = df['protrend_id_operon'].tolist()
        to_identifiers = df['genes'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToOrganismConnector(LiteratureConnector):
    default_from_node = RegulatoryInteraction
    default_to_node = Organism
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='rin',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        rin = apply_processors(rin, taxonomy=to_int_str)

        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = apply_processors(organism, ncbi_taxonomy=to_int_str)

        df = pd.merge(rin, organism, left_on='taxonomy', right_on='ncbi_taxonomy',
                      suffixes=('_rin', '_organism'))

        from_identifiers = df['protrend_id_rin'].tolist()
        to_identifiers = df['protrend_id_organism'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EffectorToOrganismConnector(LiteratureConnector):
    default_from_node = Effector
    default_to_node = Organism
    default_connect_stack = {'effector': 'integrated_effector.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        effector = read_from_stack(stack=self._connect_stack, file='effector',
                                   default_columns=EffectorTransformer.columns, reader=read_json_frame)
        effector = apply_processors(effector, taxonomy=to_int_str)

        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = apply_processors(organism, ncbi_taxonomy=to_int_str)

        df = pd.merge(effector, organism, left_on='taxonomy', right_on='ncbi_taxonomy',
                      suffixes=('_effector', '_organism'))

        from_identifiers = df['protrend_id_effector'].tolist()
        to_identifiers = df['protrend_id_organism'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
