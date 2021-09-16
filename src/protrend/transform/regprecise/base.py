from protrend.transform import Transformer, Connector


class RegPreciseTransformer(Transformer):
    default_source: str = 'regprecise'
    default_version: str = '0.0.0'


class RegPreciseConnector(Connector):
    default_source: str = 'regprecise'
    default_version: str = '0.0.0'
