from datetime import datetime

import pytz
from scrapy.exporters import JsonLinesItemExporter

from protrend.extract.items.regprecise import (TaxonomyItem, GenomeItem, TranscriptionFactorItem, RegulogItem,
                                               TranscriptionFactorFamilyItem, RNAFamilyItem, EffectorItem, PathwayItem,
                                               RegulonItem, OperonItem, GeneItem, TFBSItem, DatabaseItem)
from protrend.extract.pipelines.json_pipeline import JSONPipeline, build_json_exporters


class RegPrecisePipeline(JSONPipeline):

    def open_spider(self, spider):

        self.process_db_item()

        self.items_types = (TaxonomyItem,
                            GenomeItem,
                            TranscriptionFactorItem,
                            RegulogItem,
                            TranscriptionFactorFamilyItem,
                            RNAFamilyItem,
                            EffectorItem,
                            PathwayItem,
                            RegulonItem,
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

    def process_db_item(self):

        db_item = DatabaseItem(name='regprecise',
                               url='https://regprecise.lbl.gov/collections.jsp',
                               doi='10.1186/1471-2164-14-745',
                               authors=['Pavel S Novichkov', 'Dmitry A Rodionov'],
                               description='Collection of Manually Curated Inferences of '
                                           'Regulons in Prokaryotic Genomes',
                               version=self.version,
                               created=datetime.utcnow().replace(tzinfo=pytz.utc))

        file = open(fr'{self.staging_area}/Database.json', 'wb')
        exporter = JsonLinesItemExporter(file)
        exporter.start_exporting()
        exporter.export_item(db_item)
        exporter.finish_exporting()
        file.close()

    def process_item(self, item, spider):

        item_name = item.__class__.__name__
        item_name = item_name.replace('Item', '')

        exporter, _ = self.exporters.get(item_name, (None, None))

        exporter.export_item(item)

        return item
