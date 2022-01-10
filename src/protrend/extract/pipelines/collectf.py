from datetime import datetime

import pytz
from scrapy.exporters import JsonLinesItemExporter

from protrend.extract.items.collectf import (TaxonomyItem,
                                             OrganismItem,
                                             TranscriptionFactorItem,
                                             RegulonItem,
                                             OperonItem,
                                             GeneItem,
                                             TFBSItem,
                                             ExperimentalEvidenceItem,
                                             CollecTFItem, DatabaseItem)
from protrend.extract.pipelines.json_pipeline import JSONPipeline, build_json_exporters


class CollecTFPipeline(JSONPipeline):

    def open_spider(self, spider):

        db_item = DatabaseItem(name='collectf',
                               url='http://www.collectf.org/browse/browse/',
                               doi='10.1093/nar/gkt1123',
                               authors=['Sefa Kılıç', 'Ivan Erill'],
                               description='CollecTF: a database of experimentally validated transcription '
                                           'factor-binding sites in Bacteria',
                               version=self.version,
                               created=datetime.utcnow().replace(tzinfo=pytz.utc))

        file = open(fr'{self.staging_area}\Database.json', 'wb')
        exporter = JsonLinesItemExporter(file)
        exporter.start_exporting()
        exporter.export_item(db_item)
        exporter.finish_exporting()
        file.close()

        self.items_types = [TaxonomyItem,
                            OrganismItem,
                            TranscriptionFactorItem,
                            RegulonItem,
                            OperonItem,
                            GeneItem,
                            TFBSItem,
                            ExperimentalEvidenceItem]

        self.exporters = build_json_exporters(self.staging_area, self.items_types)

        for exporter, _ in self.exporters.values():
            exporter.start_exporting()

    def close_spider(self, spider):

        for exporter, file in self.exporters.values():
            exporter.finish_exporting()
            file.close()

    def process_item(self, item, spider):

        if isinstance(item, CollecTFItem):

            for item_container in item.values():

                for it in item_container:
                    item_name = it.__class__.__name__
                    item_name = item_name.replace('Item', '')

                    exporter, _ = self.exporters.get(item_name, (None, None))
                    exporter.export_item(it)

        else:
            item_name = item.__class__.__name__
            item_name = item_name.replace('Item', '')

            exporter, _ = self.exporters.get(item_name, (None, None))

            exporter.export_item(item)

        return item
