from typing import Dict, Tuple

from protrend.model.model import Organism, Pathway, RegulatoryFamily, Effector, Source, Publication
from protrend.transform.settings import TransformerSettings


class RegPreciseSettings(TransformerSettings):
    default_source: str = 'regprecise'


class SourceSettings(RegPreciseSettings):
    default_node: Source = Source
    default_node_factors: Tuple[str] = ('name', )


class OrganismSettings(RegPreciseSettings):
    default_files: Dict[str, str] = {'genome': 'Genome.json'}
    default_node: Organism = Organism
    default_node_factors: Tuple[str] = ('ncbi_taxonomy', 'name')


class EffectorSettings(RegPreciseSettings):
    default_files: Dict[str, str] = {'effector': 'Effector.json'}
    default_node: Effector = Effector
    default_node_factors: Tuple[str] = ('name', )


class PathwaySettings(RegPreciseSettings):
    default_files: Dict[str, str] = {'pathway': 'Pathway.json'}
    default_node: Pathway = Pathway
    default_node_factors: Tuple[str] = ('name', )


class RegulatoryFamilySettings(RegPreciseSettings):
    default_files: Dict[str, str] = {'tf_family': 'TranscriptionFactorFamily.json',
                                     'tf': 'TranscriptionFactor.json',
                                     'rna': 'RNAFamily.json'}
    default_node: RegulatoryFamily = RegulatoryFamily
    default_node_factors: Tuple[str] = ('name', 'rfam')


class PublicationSettings(RegPreciseSettings):
    default_files: Dict[str, str] = {'tf_family': 'TranscriptionFactorFamily.json',
                                     'tf': 'TranscriptionFactor.json',
                                     'rna': 'RNAFamily.json'}
    default_node: RegulatoryFamily = Publication
    default_node_factors: Tuple[str] = ('pmid',)