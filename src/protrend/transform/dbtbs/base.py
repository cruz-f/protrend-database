from protrend.transform import Transformer, Connector


class DBTBSTransformer(Transformer):
    default_source = 'dbtbs'
    default_version = '0.0.3'


class DBTBSConnector(Connector):
    default_source: str = 'dbtbs'
    default_version: str = '0.0.3'
