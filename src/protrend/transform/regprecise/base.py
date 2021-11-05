from abc import abstractmethod

from protrend.transform import Transformer, Connector


class RegPreciseTransformer(Transformer, source='regprecise', version='0.0.0', register=False):

    @abstractmethod
    def transform(self):
        pass


class RegPreciseConnector(Connector, source='regprecise', version='0.0.0', register=False):

    @abstractmethod
    def connect(self):
        pass
