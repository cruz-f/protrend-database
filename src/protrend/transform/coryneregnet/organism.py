import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import Organism, Regulator, Operon, Gene, TFBS, RegulatoryInteraction
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer, CoryneRegNetConnector
from protrend.transform.coryneregnet.gene import GeneTransformer
from protrend.transform.coryneregnet.operon import OperonTransformer
from protrend.transform.coryneregnet.regulator import RegulatorTransformer
from protrend.transform.coryneregnet.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.transform.coryneregnet.tfbs import TFBSTransformer
from protrend.transform.processors import apply_processors, to_set_list, to_int_str
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


class RegulatorToOrganismConnector(CoryneRegNetConnector):
    default_from_node = Regulator
    default_to_node = Organism
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)

        df = pd.merge(regulator, organism, left_on='taxonomy', right_on='ncbi_taxonomy',
                      suffixes=('_regulator', '_organism'))

        from_identifiers = df['protrend_id_regulator'].tolist()
        to_identifiers = df['protrend_id_organism'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class OperonToOrganismConnector(CoryneRegNetConnector):
    default_from_node = Operon
    default_to_node = Organism
    default_connect_stack = {'operon': 'integrated_operon.json', 'gene': 'integrated_gene.json',
                             'organism': 'integrated_organism.json'}

    def connect(self):
        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        operon = apply_processors(operon, genes=to_set_list)
        operon = operon.explode(column='genes')

        operon_gene = pd.merge(operon, gene, left_on='genes', right_on='protrend_id', suffixes=('_operon', '_gene'))
        operon_gene = apply_processors(operon_gene, taxonomy=to_int_str)

        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = apply_processors(organism, ncbi_taxonomy=to_int_str)

        df = pd.merge(organism, operon_gene, left_on='ncbi_taxonomy', right_on='taxonomy')

        df = OrganismTransformer.drop_duplicates(df=df, subset=['protrend_id', 'protrend_id_operon'],
                                                 perfect_match=True, preserve_nan=True)
        df = df.dropna(subset=['protrend_id', 'protrend_id_operon'])

        from_identifiers = df['protrend_id_operon'].tolist()
        to_identifiers = df['protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class GeneToOrganismConnector(CoryneRegNetConnector):
    default_from_node = Gene
    default_to_node = Organism
    default_connect_stack = {'gene': 'integrated_gene.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = apply_processors(gene, taxonomy=to_int_str)

        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = apply_processors(organism, ncbi_taxonomy=to_int_str)

        df = pd.merge(organism, gene, left_on='ncbi_taxonomy', right_on='taxonomy', suffixes=('_organism', '_gene'))

        from_identifiers = df['protrend_id_gene'].tolist()
        to_identifiers = df['protrend_id_organism'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class TFBSToOrganismConnector(CoryneRegNetConnector):
    default_from_node = TFBS
    default_to_node = Organism
    default_connect_stack = {'tfbs': 'integrated_tfbs.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        tfbs = apply_processors(tfbs, taxonomy=to_int_str)

        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = apply_processors(organism, ncbi_taxonomy=to_int_str)

        df = pd.merge(organism, tfbs, left_on='ncbi_taxonomy', right_on='taxonomy', suffixes=('_organism', '_tfbs'))

        from_identifiers = df['protrend_id_tfbs'].tolist()
        to_identifiers = df['protrend_id_organism'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToOrganismConnector(CoryneRegNetConnector):
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

        df = pd.merge(organism, rin, left_on='ncbi_taxonomy', right_on='taxonomy', suffixes=('_organism', '_rin'))

        from_identifiers = df['protrend_id_rin'].tolist()
        to_identifiers = df['protrend_id_organism'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
