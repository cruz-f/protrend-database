from protrend.transform import Transformer, Connector


class CollectfTransformer(Transformer):
    default_source: str = 'collectf'
    default_version: str = '0.0.1'


class CollectfConnector(Connector):
    default_source: str = 'collectf'
    default_version: str = '0.0.1'
