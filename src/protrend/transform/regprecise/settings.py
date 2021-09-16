from protrend.model.model import Organism, Pathway, RegulatoryFamily, Effector, Source, Publication, Regulator, Gene, \
    Operon, TFBS, RegulatoryInteraction
from protrend.transform.settings import ConnectorSettings


# ----------------------------------------------------------------------
# Connections
# ----------------------------------------------------------------------
class RegPreciseConnections(ConnectorSettings):
    default_source: str = 'regprecise'


class EffectorToSource(RegPreciseConnections):
    default_from_node = Effector
    default_to_node = Source
    default_connect = {'effector': 'integrated_effector.json',
                                       'source': 'integrated_source.json'}


class GeneToSource(RegPreciseConnections):
    default_from_node = Gene
    default_to_node = Source
    default_connect = {'gene': 'integrated_gene.json',
                                       'source': 'integrated_source.json'}


class OperonToSource(RegPreciseConnections):
    default_from_node = Operon
    default_to_node = Source
    default_connect = {'operon': 'integrated_operon.json',
                                       'source': 'integrated_source.json'}


class OrganismToSource(RegPreciseConnections):
    default_from_node = Organism
    default_to_node = Source
    default_connect = {'organism': 'integrated_organism.json',
                                       'source': 'integrated_source.json'}


class PathwayToSource(RegPreciseConnections):
    default_from_node = Pathway
    default_to_node = Source
    default_connect = {'pathway': 'integrated_pathway.json',
                                       'source': 'integrated_source.json'}


class RegulatorToSource(RegPreciseConnections):
    default_from_node = Regulator
    default_to_node = Source
    default_connect = {'regulator': 'integrated_regulator.json',
                                       'source': 'integrated_source.json'}


class RegulatoryFamilyToSource(RegPreciseConnections):
    default_from_node = RegulatoryFamily
    default_to_node = Source
    default_connect = {'regulatory_family': 'integrated_regulatoryfamily.json',
                                       'source': 'integrated_source.json'}


class TFBSToSource(RegPreciseConnections):
    default_from_node = TFBS
    default_to_node = Source
    default_connect = {'tfbs': 'integrated_tfbs.json',
                                       'source': 'integrated_source.json'}


class RegulatoryInteractionToSource(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = Source
    default_connect = {'regulatory_interaction': 'integrated_regulatoryinteraction.json',
                                       'source': 'integrated_source.json'}


class RegulatoryFamilyToPublication(RegPreciseConnections):
    default_from_node = RegulatoryFamily
    default_to_node = Publication
    default_connect = {'regulatory_family': 'integrated_regulatoryfamily.json',
                                       'publication': 'integrated_publication.json'}


class RegulatorToOrganism(RegPreciseConnections):
    default_from_node = Regulator
    default_to_node = Organism
    default_connect = {'regulator': 'integrated_regulator.json'}


class OperonToOrganism(RegPreciseConnections):
    default_from_node = Operon
    default_to_node = Organism
    default_connect = {'regulator': 'integrated_regulator.json',
                                       'operon': 'integrated_operon.json'}


class GeneToOrganism(RegPreciseConnections):
    default_from_node = Gene
    default_to_node = Organism
    default_connect = {'gene': 'integrated_gene.json'}


class TFBSToOrganism(RegPreciseConnections):
    default_from_node = TFBS
    default_to_node = Organism
    default_connect = {'tfbs': 'integrated_tfbs.json',
                                       'regulator': 'integrated_regulator.json'}


class EffectorToOrganism(RegPreciseConnections):
    default_from_node = Effector
    default_to_node = Organism
    default_connect = {'effector': 'integrated_effector.json',
                                       'regulator': 'integrated_regulator.json'}


class PathwayToRegulator(RegPreciseConnections):
    default_from_node = Pathway
    default_to_node = Regulator
    default_connect = {'pathway': 'integrated_pathway.json',
                                       'regulator': 'integrated_regulator.json'}


class PathwayToGene(RegPreciseConnections):
    default_from_node = Pathway
    default_to_node = Gene
    default_connect = {'pathway': 'integrated_pathway.json',
                                       'regulator': 'integrated_regulator.json',
                                       'gene': 'integrated_gene.json'}


class RegulatoryFamilyToRegulator(RegPreciseConnections):
    default_from_node = RegulatoryFamily
    default_to_node = Regulator
    default_connect = {'regulatory_family': 'integrated_regulatoryfamily.json',
                                       'regulator': 'integrated_regulator.json'}


class EffectorToRegulator(RegPreciseConnections):
    default_from_node = Effector
    default_to_node = Regulator
    default_connect = {'effector': 'integrated_effector.json',
                                       'regulator': 'integrated_regulator.json'}


class OperonToRegulator(RegPreciseConnections):
    default_from_node = Operon
    default_to_node = Regulator
    default_connect = {'operon': 'integrated_operon.json',
                                       'regulator': 'integrated_regulator.json'}


class OperonToGene(RegPreciseConnections):
    default_from_node = Operon
    default_to_node = Gene
    default_connect = {'operon': 'integrated_operon.json'}


class OperonToTFBS(RegPreciseConnections):
    default_from_node = Operon
    default_to_node = TFBS
    default_connect = {'operon': 'integrated_operon.json'}


class GeneToTFBS(RegPreciseConnections):
    default_from_node = Gene
    default_to_node = TFBS
    default_connect = {'operon': 'integrated_operon.json'}


class GeneToRegulator(RegPreciseConnections):
    default_from_node = Gene
    default_to_node = Regulator
    default_connect = {'operon': 'integrated_operon.json',
                                       'regulator': 'integrated_regulator.json'}


class TFBSToRegulator(RegPreciseConnections):
    default_from_node = TFBS
    default_to_node = Regulator
    default_connect = {'operon': 'integrated_operon.json',
                                       'regulator': 'integrated_regulator.json'}


class RegulatoryInteractionToOrganism(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = Organism
    default_connect = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}


class RegulatoryInteractionToEffector(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = Effector
    default_connect = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}


class RegulatoryInteractionToRegulator(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = Regulator
    default_connect = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}


class RegulatoryInteractionToOperon(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = Operon
    default_connect = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}


class RegulatoryInteractionToGene(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = Gene
    default_connect = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}


class RegulatoryInteractionToTFBS(RegPreciseConnections):
    default_from_node = RegulatoryInteraction
    default_to_node = TFBS
    default_connect = {'regulatory_interaction': 'integrated_regulatoryinteraction.json'}
