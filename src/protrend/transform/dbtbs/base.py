from abc import abstractmethod

from protrend.transform import Transformer, Connector


class DBTBSTransformer(Transformer, source='dbtbs', version='0.0.3', register=False):

    @abstractmethod
    def transform(self):
        pass


class DBTBSConnector(Connector, source='dbtbs', version='0.0.3', register=False):

    @abstractmethod
    def connect(self):
        pass
