from typing import Dict, Tuple

from protrend.model.model import Organism, Pathway, RegulatoryFamily, Effector, Source, Publication, Regulator, Gene, \
    Operon, Regulon, RegulatoryInteraction, TFBS
from protrend.transform.settings import TransformerSettings


class RegPreciseSettings(TransformerSettings):
    default_source: str = 'regprecise'


class EffectorSettings(RegPreciseSettings):
    default_node: Effector = Effector
    default_node_factors: Tuple[str] = ('name',)
    default_transform: Dict[str, str] = {'effector': 'Effector.json'}
    default_connect: Dict[str, str] = {'from': 'integrated_effector.csv',
                                       'to': 'integrated_source.csv'}
    default_order = 90


class GeneSettings(RegPreciseSettings):
    default_node: Gene = Gene
    default_node_factors: Tuple[str] = ('uniprot_accession', 'ncbi_protein', 'ncbi_gene',
                                        'genbank_accession', 'refseq_accession',
                                        'locus_tag')
    default_transform: Dict[str, str] = {'gene': 'Gene.json',
                                         'regulator': 'integrated_regulator.csv'}
    default_connect: Dict[str, str] = {'from': 'integrated_gene.csv',
                                       'to_source': 'integrated_source.csv',
                                       'to_organism': 'integrated_organism.csv',
                                       'to_regulator': 'integrated_regulator.csv'}
    default_order = 80


class OperonSettings(RegPreciseSettings):
    default_node: Operon = Operon
    default_node_factors: Tuple[str] = ('genes',)
    default_transform: Dict[str, str] = {'operon': 'Operon.json',
                                         'gene': 'integrated_gene.csv'}
    default_connect: Dict[str, str] = {'from': 'integrated_operon.csv',
                                       'to_source': 'integrated_source.csv',
                                       'to_organism': 'integrated_organism.csv',
                                       'to_regulator': 'integrated_regulator.csv',
                                       'to_gene': 'integrated_gene.csv',
                                       'to_tfbs': 'integrated_tfbs.csv'}
    default_order = 70


class OrganismSettings(RegPreciseSettings):
    default_node: Organism = Organism
    default_node_factors: Tuple[str] = ('ncbi_taxonomy', 'name')
    default_transform: Dict[str, str] = {'genome': 'Genome.json'}
    default_connect: Dict[str, str] = {'from': 'integrated_organism.csv',
                                       'to': 'integrated_source.csv'}
    default_order = 90


class PathwaySettings(RegPreciseSettings):
    default_node: Pathway = Pathway
    default_node_factors: Tuple[str] = ('name',)
    default_transform: Dict[str, str] = {'pathway': 'Pathway.json'}
    default_connect: Dict[str, str] = {'from': 'integrated_pathway.csv',
                                       'to': 'integrated_source.csv'}
    default_order = 90


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
                                         'organism': 'integrated_organism.csv'}
    default_connect: Dict[str, str] = {'from': 'integrated_regulator.csv',
                                       'to_source': 'integrated_source.csv',
                                       'to_organism': 'integrated_organism.csv',
                                       'to_effector': 'integrated_effector.csv',
                                       'to_pathway': 'integrated_pathway.csv',
                                       'to_regulatory_family': 'integrated_regulatory_family.csv',
                                       'to_publication': 'integrated_publication.csv',}
    default_order = 80


class RegulatoryInteractionSettings(RegPreciseSettings):
    default_node: RegulatoryInteraction = RegulatoryInteraction
    default_node_factors: Tuple[str] = ('regulator', 'operon', 'genes', 'tfbss', 'effectors', 'regulatory_effect')
    default_transform: Dict[str, str] = {'regulator': 'integrated_regulator.csv',
                                         'operon': 'integrated_operon.csv',
                                         'gene': 'integrated_gene.csv',
                                         'tfbs': 'integrated_tfbs.csv',
                                         'effector': 'integrated_effector.csv'}
    default_connect: Dict[str, str] = {'from': 'integrated_operon.csv',
                                       'to_source': 'integrated_source.csv',
                                       'to_organism': 'integrated_organism.csv',
                                       'to_effector': 'integrated_effector.csv',
                                       'to_regulator': 'integrated_regulator.csv',
                                       'to_operon': 'integrated_operon.csv',
                                       'to_gene': 'integrated_gene.csv',
                                       'to_tfbs': 'integrated_tfbs.csv'}
    default_order = 60


class RegulonSettings(RegPreciseSettings):
    default_node: Regulon = Regulon
    default_node_factors: Tuple[str] = ('regulator', 'operon', 'genes', 'tfbss')
    default_transform: Dict[str, str] = {'regulator': 'integrated_regulator.csv',
                                         'operon': 'integrated_operon.csv',
                                         'gene': 'integrated_gene.csv',
                                         'tfbs': 'integrated_tfbs.csv'}
    default_connect: Dict[str, str] = {'from': 'integrated_operon.csv',
                                       'to_source': 'integrated_source.csv',
                                       'to_organism': 'integrated_organism.csv',
                                       'to_regulator': 'integrated_regulator.csv',
                                       'to_operon': 'integrated_operon.csv',
                                       'to_gene': 'integrated_gene.csv',
                                       'to_tfbs': 'integrated_tfbs.csv'}
    default_order = 60


class RegulatoryFamilySettings(RegPreciseSettings):
    default_node: RegulatoryFamily = RegulatoryFamily
    default_node_factors: Tuple[str] = ('name', 'rfam')
    default_transform: Dict[str, str] = {'tf_family': 'TranscriptionFactorFamily.json',
                                         'tf': 'TranscriptionFactor.json',
                                         'rna': 'RNAFamily.json'}
    default_connect: Dict[str, str] = {'from': 'integrated_regulatory_family.csv',
                                       'to_source': 'integrated_source.csv',
                                       'to_publication': 'integrated_publication.csv'}
    default_order = 90


class SourceSettings(RegPreciseSettings):
    default_node: Source = Source
    default_node_factors: Tuple[str] = ('name',)
    default_order = 100


class TFBSSettings(RegPreciseSettings):
    default_node: TFBS = TFBS
    default_node_factors: Tuple[str] = ()
    default_transform: Dict[str, str] = {'tfbs': 'TFBS.json',
                                         'organism': 'integrated_organism.csv'}
    default_connect: Dict[str, str] = {'from': 'integrated_gene.csv',
                                       'to_source': 'integrated_source.csv',
                                       'to_organism': 'integrated_organism.csv',
                                       'to_regulator': 'integrated_regulator.csv',
                                       'to_operon': 'integrated_operon.csv',
                                       'to_gene': 'integrated_gene.csv'}
    default_order = 80
