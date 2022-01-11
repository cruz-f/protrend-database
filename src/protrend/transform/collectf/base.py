from abc import abstractmethod

from protrend.transform import Transformer, Connector


class CollecTFTransformer(Transformer, source='collectf', version='0.0.1', register=False):

    @abstractmethod
    def transform(self):
        pass


class CollecTFConnector(Connector, source='collectf', version='0.0.1', register=False):

    @abstractmethod
    def connect(self):
        pass
