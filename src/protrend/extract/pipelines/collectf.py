from protrend.extract.items.collectf import TaxonomyItem, OrganismItem, TranscriptionFactorItem, RegulonItem, \
    OperonItem, GeneItem, TFBSItem, ExperimentalEvidenceItem, CollecTFItem
from protrend.extract.pipelines.neo_pipeline import NeoPipeline
from protrend.models.collectf import (Database, CollecTFVersion, Taxonomy, Organism,
                                      TranscriptionFactor, Regulon, Operon, Gene, TFBS,
                                      ExperimentalEvidence)
from protrend.utils.node_importer import NodeImporter


class CollecTFPipeline(NeoPipeline):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.taxa = NodeImporter(node_cls=Taxonomy, path=self._import_folder)
        self.organisms = NodeImporter(node_cls=Organism, path=self._import_folder)
        self.tfs = NodeImporter(node_cls=TranscriptionFactor, path=self._import_folder)
        self.regulons = NodeImporter(node_cls=Regulon, path=self._import_folder)
        self.operons = NodeImporter(node_cls=Operon, path=self._import_folder)
        self.genes = NodeImporter(node_cls=Gene, path=self._import_folder)
        self.tfbs = NodeImporter(node_cls=TFBS, path=self._import_folder)
        self.exps = NodeImporter(node_cls=ExperimentalEvidence, path=self._import_folder)

    @property
    def importers(self):
        return [self.taxa, self.organisms, self.tfs,
                self.regulons, self.operons, self.genes, self.tfbs,
                self.exps]

    @property
    def version(self) -> CollecTFVersion:

        if self._version is None:
            versions = CollecTFVersion.nodes.order_by('-created')

            if versions:
                self._version = versions[0]

            else:
                self._version = CollecTFVersion(name='0.0.0').save()

        if isinstance(self._version, str):
            try:
                self._version = CollecTFVersion.nodes.get(name=self._version)

            except CollecTFVersion.DoesNotExist:

                self._version = CollecTFVersion(name=self._version).save()

        return self._version

    @property
    def database_node(self) -> Database:

        if self._database_node is None:

            node = Database.get_by_version(attr='name', value='collectf', version=self.version)

            if node is None:
                node = Database().save()

            self._database_node = node

        return self._database_node

    def close_spider(self, spider):

        cls_nodes = [importer.node_cls.cls_name()
                     for importer in self.importers
                     if importer.nodes]

        args = []
        for importer in self.importers:
            if importer.node_cls.cls_name() in cls_nodes:
                importer.build_imports()

                for arg in importer.args:
                    if importer.node_cls.cls_name() in arg:
                        args.append(arg)

        self.database.import_csv_data(*args)

        # self.database.connect()
        #
        # self.database_node.connect_to_version(self.version)

    def process_item(self, item, spider):

        if isinstance(item, CollecTFItem):

            for item_container in item.values():

                for it in item_container:

                    if isinstance(it, TaxonomyItem):

                        node = Taxonomy.from_item(item=it, save=False)
                        self.taxa.append(node)

                    elif isinstance(it, OrganismItem):

                        node = Organism.from_item(item=it, save=False)
                        self.organisms.append(node)

                    elif isinstance(it, TranscriptionFactorItem):

                        node = TranscriptionFactor.from_item(item=it, save=False)
                        self.tfs.append(node)

                    elif isinstance(it, RegulonItem):

                        node = Regulon.from_item(item=it, save=False)
                        self.regulons.append(node)

                    elif isinstance(it, OperonItem):

                        node = Operon.from_item(item=it, save=False)
                        self.operons.append(node)

                    elif isinstance(it, GeneItem):

                        node = Gene.from_item(item=it, save=False)
                        self.genes.append(node)

                    elif isinstance(it, TFBSItem):

                        node = TFBS.from_item(item=it, save=False)
                        self.tfbs.append(node)

                    elif isinstance(it, ExperimentalEvidenceItem):

                        node = ExperimentalEvidence.from_item(item=it, save=False)
                        self.exps.append(node)

        elif isinstance(item, TaxonomyItem):

            node = Taxonomy.from_item(item=item, save=False)
            self.taxa.append(node)

        elif isinstance(item, OrganismItem):

            node = Organism.from_item(item=item, save=False)
            self.organisms.append(node)

        elif isinstance(item, TranscriptionFactorItem):

            node = TranscriptionFactor.from_item(item=item, save=False)
            self.tfs.append(node)

        elif isinstance(item, RegulonItem):

            node = Regulon.from_item(item=item, save=False)
            self.regulons.append(node)

        elif isinstance(item, OperonItem):

            node = Operon.from_item(item=item, save=False)
            self.operons.append(node)

        elif isinstance(item, GeneItem):

            node = Gene.from_item(item=item, save=False)
            self.genes.append(node)

        elif isinstance(item, TFBSItem):

            node = TFBS.from_item(item=item, save=False)
            self.tfbs.append(node)

        elif isinstance(item, ExperimentalEvidenceItem):

            node = ExperimentalEvidence.from_item(item=item, save=False)
            self.exps.append(node)

        return item
