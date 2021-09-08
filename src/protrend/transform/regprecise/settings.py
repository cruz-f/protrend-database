from typing import Dict, Tuple

from protrend.model.model import Organism, Pathway, RegulatoryFamily, Effector, Source, Publication, Regulator, Gene, \
    Operon, TFBS, RegulatoryInteraction
from protrend.transform.settings import TransformerSettings, ConnectorSettings


# ----------------------------------------------------------------------
# Nodes
# ----------------------------------------------------------------------
class RegPreciseSettings(TransformerSettings):
    default_source: str = 'regprecise'


class EffectorSettings(RegPreciseSettings):
    default_node: Effector = Effector
    default_node_factors: Tuple[str] = ('name',)
    default_transform: Dict[str, str] = {'effector': 'Effector.json'}
    default_order = 100


class GeneSettings(RegPreciseSettings):
    default_node: Gene = Gene
    default_node_factors: Tuple[str] = ('uniprot_accession', 'ncbi_protein', 'ncbi_gene',
                                        'genbank_accession', 'refseq_accession',
                                        'locus_tag')
    default_transform: Dict[str, str] = {'gene': 'Gene.json',
                                         'regulator': 'integrated_regulator.json'}
    default_order = 80


class OperonSettings(RegPreciseSettings):
    default_node: Operon = Operon
    default_node_factors: Tuple[str] = ()
    default_transform: Dict[str, str] = {'operon': 'Operon.json',
                                         'gene': 'integrated_gene.json',
                                         'tfbs': 'integrated_tfbs.json'}
    default_order = 60


class OrganismSettings(RegPreciseSettings):
    default_node: Organism = Organism
    default_node_factors: Tuple[str] = ('ncbi_taxonomy', 'name')
    default_transform: Dict[str, str] = {'genome': 'Genome.json'}
    default_order = 100


class PathwaySettings(RegPreciseSettings):
    default_node: Pathway = Pathway
    default_node_factors: Tuple[str] = ('name',)
    default_transform: Dict[str, str] = {'pathway': 'Pathway.json'}
    default_order = 100


class PublicationSettings(RegPreciseSettings):
    default_node: Publication = Publication
    default_node_factors: Tuple[str] = ('pmid',)
    default_transform: Dict[str, str] = {'tf_family': 'TranscriptionFactorFamily.json',
                                         'tf': 'TranscriptionFactor.json',
                                         'rna': 'RNAFamily.json'}
    default_order = 100


class RegulatorSettings(RegPreciseSettings):
    default_node: Regulator = Regulator
    default_node_factors: Tuple[str] = ('uniprot_accession', 'ncbi_protein', 'ncbi_gene',
                                        'genbank_accession', 'refseq_accession',
                                        'locus_tag')
    default_transform: Dict[str, str] = {'regulon': 'Regulon.json',
                                         'organism': 'integrated_organism.json'}
    default_order = 90


class RegulatoryFamilySettings(RegPreciseSettings):
    default_node: RegulatoryFamily = RegulatoryFamily
    default_node_factors: Tuple[str] = ('name', 'rfam')
    default_transform: Dict[str, str] = {'tf_family': 'TranscriptionFactorFamily.json',
                                         'tf': 'TranscriptionFactor.json',
                                         'rna': 'RNAFamily.json'}
    default_order = 100


class RegulatoryInteractionSettings(RegPreciseSettings):
    default_node: RegulatoryFamily = RegulatoryInteraction
    default_node_factors: Tuple[str] = ()
    default_transform: Dict[str, str] = {'effector': 'integrated_effector.json',
                                         'regulator': 'integrated_regulator.json',
                                         'operon': 'integrated_operon.json',
                                         'gene': 'integrated_gene.json',
                                         'tfbs': 'integrated_tfbs.json'}
    default_order = 50


class SourceSettings(RegPreciseSettings):
    default_node: Source = Source
    default_node_factors: Tuple[str] = ('name',)
    default_order = 100


class TFBSSettings(RegPreciseSettings):
    default_node: TFBS = TFBS
    default_node_factors: Tuple[str] = ()
    default_transform: Dict[str, str] = {'tfbs': 'TFBS.json',
                                         'gene': 'integrated_gene.json'}
    default_order = 70


# ----------------------------------------------------------------------
# Connections
# ----------------------------------------------------------------------
class RegPreciseConnections(ConnectorSettings):
    default_source: str = 'regprecise'


class EffectorToSource(RegPreciseConnections):
    default_from_node = Effector
    default_to_node = Source
    default_connect: Dict[str, str] = {'effector': 'integrated_effector.json',
                                       'source': 'integrated_source.json'}


class GeneToSource(RegPreciseConnections):
    default_from_node = Gene
    default_to_node = Source
    default_connect: Dict[str, str] = {'gene': 'integrated_gene.json',
                                       'source': 'integrated_source.json'}


class OperonToSource(RegPreciseConnections):
    default_from_node = Operon
    default_to_node = Source
    default_connect: Dict[str, str] = {'operon': 'integrated_operon.json',
                                       'source': 'integrated_source.json'}


class OrganismToSource(RegPreciseConnections):
    default_from_node = Organism
    default_to_node = Source
    default_connect: Dict[str, str] = {'organism': 'integrated_organism.json',
                                       'source': 'integrated_source.json'}


class PathwayToSource(RegPreciseConnections):
    default_from_node = Pathway
    default_to_node = Source
    default_connect: Dict[str, str] = {'pathway': 'integrated_pathway.json',
                                       'source': 'integrated_source.json'}


class RegulatorToSource(RegPreciseConnections):
    default_from_node = Regulator
    default_to_node = Source
    default_connect: Dict[str, str] = {'regulator': 'integrated_regulator.json',
                                       'source': 'integrated_source.json'}


class RegulatoryFamilyToSource(RegPreciseConnections):
    default_from_node = RegulatoryFamily
    default_to_node = Source
    default_connect: Dict[str, str] = {'regulatory_family': 'integrated_regulatoryfamily.json',
                                       'source': 'integrated_source.json'}


class TFBSToSource(RegPreciseConnections):
    default_from_node = TFBS
    default_to_node = Source
    default_connect: Dict[str, str] = {'tfbs': 'integrated_tfbs.json',
                                       'source': 'integrated_source.json'}


class RegulatoryInteractionToSource(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = Source
    default_connect: Dict[str, str] = {'regulatory_interaction': 'integrated_regulatoryinteraction.json',
                                       'source': 'integrated_source.json'}


class RegulatoryFamilyToPublication(RegPreciseConnections):
    default_from_node = RegulatoryFamily
    default_to_node = Publication
    default_connect: Dict[str, str] = {'regulatory_family': 'integrated_regulatoryfamily.json',
                                       'publication': 'integrated_publication.json'}


class RegulatorToOrganism(RegPreciseConnections):
    default_from_node = Regulator
    default_to_node = Organism
    default_connect: Dict[str, str] = {'regulator': 'integrated_regulator.json'}


class OperonToOrganism(RegPreciseConnections):
    default_from_node = Operon
    default_to_node = Organism
    default_connect: Dict[str, str] = {'regulator': 'integrated_regulator.json',
                                       'operon': 'integrated_operon.json'}


class GeneToOrganism(RegPreciseConnections):
    default_from_node = Gene
    default_to_node = Organism
    default_connect: Dict[str, str] = {'gene': 'integrated_gene.json'}


class TFBSToOrganism(RegPreciseConnections):
    default_from_node = TFBS
    default_to_node = Organism
    default_connect: Dict[str, str] = {'tfbs': 'integrated_tfbs.json',
                                       'regulator': 'integrated_regulator.json'}


class EffectorToOrganism(RegPreciseConnections):
    default_from_node = Effector
    default_to_node = Organism
    default_connect: Dict[str, str] = {'effector': 'integrated_effector.json',
                                       'regulator': 'integrated_regulator.json'}


class PathwayToRegulator(RegPreciseConnections):
    default_from_node = Pathway
    default_to_node = Regulator
    default_connect: Dict[str, str] = {'pathway': 'integrated_pathway.json',
                                       'regulator': 'integrated_regulator.json'}


class PathwayToGene(RegPreciseConnections):
    default_from_node = Pathway
    default_to_node = Gene
    default_connect: Dict[str, str] = {'pathway': 'integrated_pathway.json',
                                       'regulator': 'integrated_regulator.json',
                                       'gene': 'integrated_gene.json'}


class RegulatoryFamilyToRegulator(RegPreciseConnections):
    default_from_node = RegulatoryFamily
    default_to_node = Regulator
    default_connect: Dict[str, str] = {'regulatory_family': 'integrated_regulatoryfamily.json',
                                       'regulator': 'integrated_regulator.json'}


class EffectorToRegulator(RegPreciseConnections):
    default_from_node = Effector
    default_to_node = Regulator
    default_connect: Dict[str, str] = {'effector': 'integrated_effector.json',
                                       'regulator': 'integrated_regulator.json'}


class OperonToRegulator(RegPreciseConnections):
    default_from_node = Operon
    default_to_node = Regulator
    default_connect: Dict[str, str] = {'operon': 'integrated_operon.json',
                                       'regulator': 'integrated_regulator.json'}


class OperonToGene(RegPreciseConnections):
    default_from_node = Operon
    default_to_node = Gene
    default_connect: Dict[str, str] = {'operon': 'integrated_operon.json'}


class OperonToTFBS(RegPreciseConnections):
    default_from_node = Operon
    default_to_node = TFBS
    default_connect: Dict[str, str] = {'operon': 'integrated_operon.json'}


class GeneToTFBS(RegPreciseConnections):
    default_from_node = Gene
    default_to_node = TFBS
    default_connect: Dict[str, str] = {'operon': 'integrated_operon.json'}


class GeneToRegulator(RegPreciseConnections):
    default_from_node = Gene
    default_to_node = Regulator
    default_connect: Dict[str, str] = {'operon': 'integrated_operon.json',
                                       'regulator': 'integrated_regulator.json'}


class TFBSToRegulator(RegPreciseConnections):
    default_from_node = TFBS
    default_to_node = Regulator
    default_connect: Dict[str, str] = {'operon': 'integrated_operon.json',
                                       'regulator': 'integrated_regulator.json'}


class RegulatoryInteractionToOrganism(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = Organism
    default_connect: Dict[str, str] = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}


class RegulatoryInteractionToEffector(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = Effector
    default_connect: Dict[str, str] = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}


class RegulatoryInteractionToRegulator(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = Regulator
    default_connect: Dict[str, str] = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}


class RegulatoryInteractionToOperon(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = Operon
    default_connect: Dict[str, str] = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}


class RegulatoryInteractionToGene(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = Gene
    default_connect: Dict[str, str] = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}


class RegulatoryInteractionToTFBS(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = TFBS
    default_connect: Dict[str, str] = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}
