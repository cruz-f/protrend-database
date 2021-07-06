from typing import List, Type, Dict, Tuple, TextIO

from scrapy import Item
from scrapy.exporters import JsonLinesItemExporter


class JSONPipeline:

    def __init__(self,
                 staging_area: str,
                 version: str):
        self.staging_area: str = staging_area
        self.version: str = version

        self.exporters: Dict[str, Tuple[JsonLinesItemExporter, TextIO]] = {}
        self.items_types: Tuple[Type[Item]] = ()

    @classmethod
    def from_crawler(cls, crawler):
        staging_area = crawler.settings.get('staging_area')
        version = crawler.settings.get('version')

        return cls(staging_area=staging_area,
                   version=version)

    def open_spider(self, spider):
        pass

    def close_spider(self, spider):
        pass

    def process_item(self, item, spider):
        return item


def build_json_exporters(path: str, items_types: List[Type[Item]]) -> Dict[str, Tuple[JsonLinesItemExporter, TextIO]]:
    exporters = {}

    for item in items_types:
        item_name = item.__name__
        item_name = item_name.replace('Item', '')

        file = open(fr'{path}\{item_name}.json', 'wb')
        exporter = JsonLinesItemExporter(file)
        exporters[item_name] = (exporter, file)

    return exporters
