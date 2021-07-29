from typing import Dict


class RegPreciseTransformSettings:
    source: str = 'regprecise'
    version: str = '0.0.0'
    organism: Dict[str, str] = {'genome': 'Genome.json'}
    pathway: Dict[str, str] = {'pathway': 'Pathway.json'}
    effector: Dict[str, str] = {'effector': 'Effector.json'}
    regulatory_family: Dict[str, str] = {'tf_family': 'TranscriptionFactorFamily.json',
                                         'tf': 'TranscriptionFactor.json',
                                         'rna': 'RNAFamily.json'}
