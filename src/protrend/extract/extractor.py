import sys
from abc import ABCMeta, abstractmethod

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
    default_spider: str = ''
    default_version: str = '0.0.0'

    def __init__(self,
                 spider: str = None,
                 version: str = None):
        """
        The extractor object calls the respective scrapy spider for a given data source.

        :type spider: str
        :type version: str

        :param spider: The name of the data source and respective spider (e.g. regprecise, collectf, etc)
        :param version: The version of the data source in the staging area (e.g. 0.0.0, 0.0.1, etc)
        """

        self._spider = spider
        self._version = version

    # --------------------------------------------------------
    # Static properties
    # --------------------------------------------------------
    @property
    def spider(self) -> str:
        if not self._spider:
            return self.default_spider

        return self._spider

    @property
    def version(self) -> str:
        if not self._version:
            return self.default_version

        return self._version

    def extract(self):
        src_path = Settings.ROOT_PATH.parent
        sys.path.insert(0, str(src_path))
        from .run_spider import run_spider
        return run_spider(spider=self.spider, staging_area=Settings.STAGING_AREA_PATH, version=self.version)


class RegPreciseExtractor(Extractor):
    default_spider = 'regprecise'
    default_version = '0.0.0'


class CollectfExtractor(Extractor):
    default_spider = 'collectf'
    default_version = '0.0.1'


class DBTBSExtractor(Extractor):
    default_spider = 'dbtbs'
    default_version = '0.0.3'
