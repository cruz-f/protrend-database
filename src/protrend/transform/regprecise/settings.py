from typing import Dict


class RegPreciseTransformSettings:
    source: str = 'regprecise'
    version: str = '0.0.0'
    organism: Dict[str, str] = {'taxonomy': 'Taxonomy.json',
                                'genome': 'Genome.json'}
