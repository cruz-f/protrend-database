from protrend.transform import Transformer, Connector
from protrend.model import Source, Motif
from protrend.transform.mix_ins import SourceMixIn
from protrend.utils import SetList
from protrend.utils.constants import CURATION


class SourceTransformer(SourceMixIn, Transformer,
                        source='functional_tfbs',
                        version='0.0.0',
                        node=Source,
                        order=100,
                        register=True):
    name = ['functional_tfbs']
    type = [CURATION]
    url = ['https://protrend.bio.di.uminho.pt/']
    doi = ['10.1016/j.csbj.2020.05.015']
    authors = [['ProTReND Team']]
    description = ['Prokaryotic Transcriptional Regulatory Networks Database']

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])


class SourceToMotifConnector(Connector,
                             source='functional_tfbs',
                             version='0.0.0',
                             from_node=Source,
                             to_node=Motif,
                             register=True):

    def connect(self):
        df = self.create_connection(source='source', target='motif', cardinality='one_to_many')
        self.stack_connections(df)
