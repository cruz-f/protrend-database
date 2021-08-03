from typing import Dict, Tuple

from protrend.model.model import Organism, Pathway, RegulatoryFamily, Effector, Source, Publication
from protrend.transform.settings import TransformerSettings, RelationshipSettings


class RegPreciseSettings(TransformerSettings):
    default_source: str = 'regprecise'


class SourceSettings(RegPreciseSettings):
    default_node: Source = Source
    default_node_factors: Tuple[str] = ('name', )


class OrganismSettings(RegPreciseSettings):
    default_files: Dict[str, str] = {'genome': 'Genome.json'}
    default_node: Organism = Organism
    default_node_factors: Tuple[str] = ('ncbi_taxonomy', 'name')
    default_relationships = [RelationshipSettings(from_node=Organism,
                                                  to_node=Source,
                                                  from_property=Organism.identifying_property,
                                                  to_property='name',
                                                  **{Organism.identifying_property: 'from_identifier',
                                                     'source_db': 'to_identifier',
                                                     'url': 'url',
                                                     'genome_id': 'external_identifier',
                                                     'api_key': 'name'})]


class EffectorSettings(RegPreciseSettings):
    default_files: Dict[str, str] = {'effector': 'Effector.json'}
    default_node: Effector = Effector
    default_node_factors: Tuple[str] = ('name', )
    default_relationships = [RelationshipSettings(from_node=Effector,
                                                  to_node=Source,
                                                  from_property=Effector.identifying_property,
                                                  to_property='name',
                                                  **{Effector.identifying_property: 'from_identifier',
                                                     'source_db': 'to_identifier',
                                                     'url': 'url',
                                                     'genome_id': 'external_identifier',
                                                     'api_key': 'name'})]


class PathwaySettings(RegPreciseSettings):
    default_files: Dict[str, str] = {'pathway': 'Pathway.json'}
    default_node: Pathway = Pathway
    default_node_factors: Tuple[str] = ('name', )
    default_relationships = [RelationshipSettings(from_node=Pathway,
                                                  to_node=Source,
                                                  from_property=Pathway.identifying_property,
                                                  to_property='name',
                                                  **{Pathway.identifying_property: 'from_identifier',
                                                     'source_db': 'to_identifier',
                                                     'url': 'url',
                                                     'genome_id': 'external_identifier',
                                                     'api_key': 'name'})]


class RegulatoryFamilySettings(RegPreciseSettings):
    default_files: Dict[str, str] = {'tf_family': 'TranscriptionFactorFamily.json',
                                     'tf': 'TranscriptionFactor.json',
                                     'rna': 'RNAFamily.json'}
    default_node: RegulatoryFamily = RegulatoryFamily
    default_node_factors: Tuple[str] = ('name', 'rfam')
    default_relationships = [RelationshipSettings(from_node=RegulatoryFamily,
                                                  to_node=Source,
                                                  from_property=RegulatoryFamily.identifying_property,
                                                  to_property='name',
                                                  **{RegulatoryFamily.identifying_property: 'from_identifier',
                                                     'source_db': 'to_identifier',
                                                     'url': 'url',
                                                     'external_identifier': 'external_identifier',
                                                     'api_key': 'name'}),
                             RelationshipSettings(from_node=RegulatoryFamily,
                                                  to_node=Publication,
                                                  from_property=RegulatoryFamily.identifying_property,
                                                  to_property='pmid',
                                                  **{RegulatoryFamily.identifying_property: 'from_identifier',
                                                     'pubmed': 'to_identifier'}),
                             ]
