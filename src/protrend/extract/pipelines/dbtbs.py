from datetime import datetime

import pytz
from scrapy.exporters import JsonLinesItemExporter

from protrend.extract.items.dbtbs import (DatabaseItem,
                                          TranscriptionFactorItem,
                                          OperonItem,
                                          GeneItem,
                                          TFBSItem)
from protrend.extract.pipelines.json_pipeline import JSONPipeline, build_json_exporters


class DBTBSPipeline(JSONPipeline):

    def open_spider(self, spider):

        db_item = DatabaseItem(name='dbtbs',
                               url='https://dbtbs.hgc.jp/',
                               doi='10.1093/nar/gkm910',
                               authors=['Nicolas Sierro', 'Kenta Nakai'],
                               description='DBTBS: a database of transcriptional regulation in Bacillus subtilis '
                                           'containing upstream intergenic conservation information',
                               version=self.version,
                               created=datetime.utcnow().replace(tzinfo=pytz.utc))

        file = open(fr'{self.staging_area}\Database.json', 'wb')
        exporter = JsonLinesItemExporter(file)
        exporter.start_exporting()
        exporter.export_item(db_item)
        exporter.finish_exporting()
        file.close()

        self.items_types = (DatabaseItem,
                            TranscriptionFactorItem,
                            OperonItem,
                            GeneItem,
                            TFBSItem)

        self.exporters = build_json_exporters(self.staging_area, self.items_types)

        for exporter, _ in self.exporters.values():
            exporter.start_exporting()

    def close_spider(self, spider):

        for exporter, file in self.exporters.values():
            exporter.finish_exporting()
            file.close()

    def process_item(self, item, spider):

        item_name = item.__class__.__name__
        item_name = item_name.replace('Item', '')

        exporter, _ = self.exporters.get(item_name, (None, None))

        exporter.export_item(item)

        return item
