from abc import abstractmethod

from protrend.transform import Transformer, Connector


class CollectfTransformer(Transformer, source='collectf', version='0.0.1', register=False):

    @abstractmethod
    def transform(self):
        pass


class CollectfConnector(Connector, source='collectf', version='0.0.1', register=False):

    @abstractmethod
    def connect(self):
        pass
