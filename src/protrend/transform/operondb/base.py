from abc import abstractmethod

from protrend.transform import Transformer, Connector


class OperonDBTransformer(Transformer, source='operondb', version='0.0.0', register=False):

    @abstractmethod
    def transform(self):
        pass


class OperonDBConnector(Connector, source='operondb', version='0.0.0', register=False):

    @abstractmethod
    def connect(self):
        pass
