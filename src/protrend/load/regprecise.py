from typing import List

from protrend.load.settings import LoaderSettings
from protrend.model.model import Organism


class RegPreciseSettings(LoaderSettings):
    default_source: str = 'regprecise'


class OrganismSettings(RegPreciseSettings):
    default_files: List[str] = ['integrated_organism.csv', 'connected_organism_source.csv']
    default_node: Organism = Organism