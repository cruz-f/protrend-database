import sys
from abc import ABCMeta, abstractmethod

from protrend.utils import Settings
from protrend.utils import DefaultProperty


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
    source = DefaultProperty('')
    version = DefaultProperty('')

    def __init_subclass__(cls, **kwargs):

        source = kwargs.get('source')
        cls.source.set_default(source)

        version = kwargs.get('version')
        cls.version.set_default(source)

        register = kwargs.pop('register', False)

        if register:
            from protrend.pipeline import Pipeline
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
        src_path = Settings.ROOT_PATH.parent
        sys.path.insert(0, str(src_path))
        from .run_spider import run_spider
        return run_spider(spider=self.source, staging_area=Settings.STAGING_AREA_PATH, version=self.version)
