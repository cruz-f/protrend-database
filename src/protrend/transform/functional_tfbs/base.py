from abc import abstractmethod

from protrend.transform import Transformer, Connector


class FunctionalTFBSTransformer(Transformer, register=False):

    @abstractmethod
    def transform(self):
        pass


class FunctionalTFBSConnector(Connector, register=False):

    @abstractmethod
    def connect(self):
        pass
