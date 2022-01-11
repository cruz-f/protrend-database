from abc import abstractmethod
from functools import partial
from typing import List, Union

from protrend.io import read_txt
from protrend.transform import Transformer, Connector
from protrend.utils import SetList


def regulondb_reader(skiprows: int, names: Union[SetList, List[str]]) -> partial:
    return partial(read_txt, skiprows=skiprows, names=names)


class RegulonDBTransformer(Transformer, source='regulondb', version='0.0.0', register=False):

    @abstractmethod
    def transform(self):
        pass


class RegulonDBConnector(Connector, source='regulondb', version='0.0.0', register=False):

    @abstractmethod
    def connect(self):
        pass
