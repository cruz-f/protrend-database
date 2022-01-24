from typing import List, Type, Dict, Tuple, TextIO, BinaryIO

from scrapy import Item
from scrapy.exporters import JsonLinesItemExporter


class JSONPipeline:

    def __init__(self,
                 data_lake: str,
                 version: str):
        self.data_lake: str = data_lake
        self.version: str = version

        self.exporters: Dict[str, Tuple[JsonLinesItemExporter, TextIO]] = {}
        self.items_types: Tuple = ()

    @classmethod
    def from_crawler(cls, crawler):
        data_lake = crawler.settings.get('data_lake')
        version = crawler.settings.get('version')

        return cls(data_lake=data_lake,
                   version=version)

    def open_spider(self, spider):
        pass

    def close_spider(self, spider):
        pass

    def process_item(self, item, spider):
        return item


def build_json_exporters(path: str, items_types: List[Type[Item]]) -> Dict[str, Tuple[JsonLinesItemExporter, BinaryIO]]:
    exporters = {}

    for item in items_types:
        item_name = item.__name__
        item_name = item_name.replace('Item', '')

        file = open(fr'{path}\{item_name}.json', 'wb')
        exporter = JsonLinesItemExporter(file)
        exporters[item_name] = (exporter, file)

    return exporters
