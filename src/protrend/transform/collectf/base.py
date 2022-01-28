from abc import abstractmethod

from protrend.transform import Transformer, Connector


class CollecTFTransformer(Transformer, register=False):

    @abstractmethod
    def transform(self):
        pass


class CollecTFConnector(Connector, register=False):

    @abstractmethod
    def connect(self):
        pass
