from protrend.model import Source
from protrend.transform.standardizer.base import StandardizerTransformer
from protrend.transform.mix_ins import SourceMixIn
from protrend.utils import SetList
from protrend.utils.constants import CURATION


class SourceTransformer(SourceMixIn, StandardizerTransformer,
                        source='standardizer',
                        version='0.0.0',
                        node=Source,
                        order=90,
                        register=True):
    name = ['standardizer']
    type = [CURATION]
    url = ['https://protrend.bio.di.uminho.pt/']
    doi = ['10.1016/j.csbj.2020.05.015']
    authors = [['ProTReND Team']]
    description = ['Prokaryotic Transcriptional Regulatory Networks Database']

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])