from abc import ABCMeta, abstractmethod


class Transformer(metaclass=ABCMeta):

    """
    Transformer interface.

    The following methods must be implemented to set up a transformer for each node
    and relationship

    A transformer is responsible for reading, processing, integrating and writing node and relationships files.

    A transformer starts with data from the staging area and ends with structured nodes and relationships.
    """

    @abstractmethod
    def read(self):
        pass

    @abstractmethod
    def process(self):
        pass

    @abstractmethod
    def integrate(self):
        pass

    @abstractmethod
    def write(self):
        pass
