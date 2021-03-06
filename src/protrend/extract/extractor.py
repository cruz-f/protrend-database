import sys
from abc import ABCMeta, abstractmethod

from protrend.extract.run_spider import run_spider
from protrend.utils import DefaultProperty
from protrend.utils import Settings


class AbstractExtractor(metaclass=ABCMeta):
    """
    Extractor interface.
    The following methods must be implemented to set up a extractor for each data source
    """

    # --------------------------------------------------------
    # Extractor API
    # --------------------------------------------------------
    @abstractmethod
    def extract(self):
        """
        The method responsible for extracting data sources.

        Interface implementation with unknown signature.
        Concrete implementations are available at the extractor.

        :return:
        """

        pass


class Extractor(AbstractExtractor):
    """
    A Extractor is responsible for extracting data sources.
    """
    source = DefaultProperty()
    version = DefaultProperty()

    def __init_subclass__(cls, **kwargs):

        source = kwargs.get('source')
        cls.source.set_default(cls, source)

        version = kwargs.get('version')
        cls.version.set_default(cls, version)

        register = kwargs.pop('register', False)

        if register:
            from protrend.pipeline.pipeline import Pipeline
            Pipeline.register_extractor(cls, **kwargs)

    def __init__(self,
                 source: str = None,
                 version: str = None):
        """
        The extractor object calls the respective scrapy spider for a given data source.

        :type source: str
        :type version: str

        :param source: The name of the data source and respective spider (e.g. regprecise, collectf, etc)
        :param version: The version of the data source in the staging area (e.g. 0.0.0, 0.0.1, etc)
        """

        self.source = source
        self.version = version

    def extract(self):
        src_path = Settings.source.parent
        sys.path.insert(0, str(src_path))
        data_lake = str(Settings.data_lake)
        return run_spider(spider=self.source, data_lake=data_lake, version=self.version)
