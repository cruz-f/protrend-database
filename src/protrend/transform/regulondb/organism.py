import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model import Regulator, Organism, Operon, Gene, TFBS, Effector, RegulatoryInteraction
from protrend.transform.regulondb.base import RegulondbTransformer, RegulondbConnector
from protrend.transform.regulondb.gene import GeneTransformer
from protrend.transform.regulondb.operon import OperonTransformer
from protrend.transform.regulondb.regulator import RegulatorTransformer
from protrend.transform.regulondb.regulatory_interaction import RegulatoryInteractionTransformer
from protrend.transform.regulondb.tfbs import TFBSTransformer
from protrend.utils import SetList


class OrganismTransformer(RegulondbTransformer,
                          source='regulondb',
                          version='0.0.0',
                          node=Organism,
                          order=100,
                          register=True):
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

        self.stack_transformed_nodes(df)

        return df


class RegulatorToOrganismConnector(RegulondbConnector,
                                   source='regulondb',
                                   version='0.0.0',
                                   from_node=Regulator,
                                   to_node=Organism,
                                   register=True):
    default_connect_stack = {'regulator': 'integrated_regulator.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)

        from_identifiers = regulator['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = organism['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class OperonToOrganismConnector(RegulondbConnector,
                                source='regulondb',
                                version='0.0.0',
                                from_node=Operon,
                                to_node=Organism,
                                register=True):
    default_connect_stack = {'operon': 'integrated_operon.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        operon = read_from_stack(stack=self._connect_stack, file='operon',
                                 default_columns=OperonTransformer.columns, reader=read_json_frame)
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)

        from_identifiers = operon['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = organism['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class GeneToOrganismConnector(RegulondbConnector,
                              source='regulondb',
                              version='0.0.0',
                              from_node=Gene,
                              to_node=Organism,
                              register=True):
    default_connect_stack = {'gene': 'integrated_gene.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)

        from_identifiers = gene['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = organism['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class TFBSToOrganismConnector(RegulondbConnector,
                              source='regulondb',
                              version='0.0.0',
                              from_node=TFBS,
                              to_node=Organism,
                              register=True):
    default_connect_stack = {'tfbs': 'integrated_tfbs.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        tfbs = read_from_stack(stack=self._connect_stack, file='tfbs',
                               default_columns=TFBSTransformer.columns, reader=read_json_frame)
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)

        from_identifiers = tfbs['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = organism['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class EffectorToOrganismConnector(RegulondbConnector,
                                  source='regulondb',
                                  version='0.0.0',
                                  from_node=Effector,
                                  to_node=Organism,
                                  register=True):
    default_connect_stack = {'effector': 'integrated_effector.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        effector = read_from_stack(stack=self._connect_stack, file='effector',
                                   default_columns=TFBSTransformer.columns, reader=read_json_frame)
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)

        from_identifiers = effector['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = organism['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)


class RegulatoryInteractionToOrganismConnector(RegulondbConnector,
                                               source='regulondb',
                                               version='0.0.0',
                                               from_node=RegulatoryInteraction,
                                               to_node=Organism,
                                               register=True):
    default_connect_stack = {'regulatory_interaction': 'integrated_regulatoryinteraction.json',
                             'organism': 'integrated_organism.json'}

    def connect(self):
        rin = read_from_stack(stack=self._connect_stack, file='regulatory_interaction',
                              default_columns=RegulatoryInteractionTransformer.columns, reader=read_json_frame)
        organism = read_from_stack(stack=self._connect_stack, file='organism',
                                   default_columns=OrganismTransformer.columns, reader=read_json_frame)

        from_identifiers = rin['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = organism['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
