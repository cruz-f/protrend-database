from abc import abstractmethod

from protrend.transform import Transformer, Connector


class OperonDBTransformer(Transformer, register=False):

    @abstractmethod
    def transform(self):
        pass


class OperonDBConnector(Connector, register=False):

    @abstractmethod
    def connect(self):
        pass
