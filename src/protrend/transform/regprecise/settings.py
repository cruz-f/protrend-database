from typing import Dict, Tuple

from protrend.model.model import Organism, Pathway, RegulatoryFamily, Effector, Source
from protrend.transform.settings import TransformerSettings, RelationshipSettings


class RegPreciseSettings(TransformerSettings):
    default_source: str = 'regprecise'


class OrganismSettings(RegPreciseSettings):
    default_files: Dict[str, str] = {'genome': 'Genome.json'}
    default_node: Organism = Organism
    default_node_factors: Tuple[str] = ('ncbi_taxonomy', 'name')
    default_relationships = [RelationshipSettings(from_node=Organism,
                                                  to_node=Source,
                                                  from_property=Organism.identifying_property,
                                                  to_property='name',
                                                  **{'url': 'url',
                                                     'genome_id': 'external_identifier',
                                                     'api_key': 'name'})]


class PathwaySettings(RegPreciseSettings):
    default_files: Dict[str, str] = {'pathway': 'Pathway.json'}
    default_node: Pathway = Pathway
    default_node_factors: Tuple[str] = ('name',)


class EffectorSettings(RegPreciseSettings):
    default_files: Dict[str, str] = {'effector': 'Effector.json'}
    default_node: Effector = Effector
    default_node_factors: Tuple[str] = ('name',)


class RegulatoryFamilySettings(RegPreciseSettings):
    default_files: Dict[str, str] = {'tf_family': 'TranscriptionFactorFamily.json',
                                     'tf': 'TranscriptionFactor.json',
                                     'rna': 'RNAFamily.json'}
    default_node: RegulatoryFamily = RegulatoryFamily
    default_node_factors: Tuple[str] = ('name',)
