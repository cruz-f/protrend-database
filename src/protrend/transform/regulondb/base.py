from protrend.transform import Transformer, Connector


class RegulondbTransformer(Transformer):
    default_source: str = 'regulondb'
    default_version: str = '0.0.0'


class RegulondbConnector(Connector):
    default_source: str = 'regulondb'
    default_version: str = '0.0.0'
