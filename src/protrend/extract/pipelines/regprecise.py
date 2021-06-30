from protrend.extract.items.regprecise import (TaxonomyItem, GenomeItem, TranscriptionFactorItem, RegulogItem,
                                               TranscriptionFactorFamilyItem, RNAFamilyItem, EffectorItem, PathwayItem,
                                               RegulonItem, OperonItem, GeneItem, TFBSItem)
from protrend.extract.pipelines.neo_pipeline import NeoPipeline
from protrend.models.regprecise import (Database, RegPreciseVersion, Taxonomy, TranscriptionFactor, TranscriptionFactorFamily,
                                        RNAFamily, Effector, Pathway, Genome, Regulog, Regulon, Operon, Gene, TFBS)
from protrend.utils.node_importer import NodeImporter


class RegPrecisePipeline(NeoPipeline):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.taxa = NodeImporter(node_cls=Taxonomy, path=self._import_folder)
        self.genomes = NodeImporter(node_cls=Genome, path=self._import_folder)
        self.tfs = NodeImporter(node_cls=TranscriptionFactor, path=self._import_folder)
        self.regulogs = NodeImporter(node_cls=Regulog, path=self._import_folder)
        self.tffams = NodeImporter(node_cls=TranscriptionFactorFamily, path=self._import_folder)
        self.rnafams = NodeImporter(node_cls=RNAFamily, path=self._import_folder)
        self.effectors = NodeImporter(node_cls=Effector, path=self._import_folder)
        self.pathways = NodeImporter(node_cls=Pathway, path=self._import_folder)
        self.regulons = NodeImporter(node_cls=Regulon, path=self._import_folder)
        self.operons = NodeImporter(node_cls=Operon, path=self._import_folder)
        self.genes = NodeImporter(node_cls=Gene, path=self._import_folder)
        self.tfbs = NodeImporter(node_cls=TFBS, path=self._import_folder)

    @property
    def importers(self):
        return [self.taxa, self.genomes, self.tfs, self.regulogs, self.tffams, self.rnafams, self.effectors,
                self.pathways, self.regulons, self.operons, self.genes, self.tfbs]

    @property
    def version(self) -> RegPreciseVersion:

        if self._version is None:
            versions = RegPreciseVersion.nodes.order_by('-created')

            if versions:
                self._version = versions[0]

            else:
                self._version = RegPreciseVersion(name='0.0.0').save()

        if isinstance(self._version, str):
            try:
                self._version = RegPreciseVersion.nodes.get(name=self._version)

            except RegPreciseVersion.DoesNotExist:

                self._version = RegPreciseVersion(name=self._version).save()

        return self._version

    @property
    def database_node(self) -> Database:

        if self._database_node is None:

            node = Database.get_by_version(attr='name', value='regprecise', version=self.version)

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

        if isinstance(item, TaxonomyItem):

            node = Taxonomy.from_item(item=item, save=False)
            self.taxa.append(node)

        elif isinstance(item, GenomeItem):

            node = Genome.from_item(item=item, save=False)
            self.genomes.append(node)

        elif isinstance(item, TranscriptionFactorItem):

            node = TranscriptionFactor.from_item(item=item, save=False)
            self.tfs.append(node)

        elif isinstance(item, RegulogItem):

            node = Regulog.from_item(item=item, save=False)
            self.regulogs.append(node)

        elif isinstance(item, TranscriptionFactorFamilyItem):

            node = TranscriptionFactorFamily.from_item(item=item, save=False)
            self.tffams.append(node)

        elif isinstance(item, RNAFamilyItem):

            node = RNAFamily.from_item(item=item, save=False)
            self.rnafams.append(node)

        elif isinstance(item, EffectorItem):

            node = Effector.from_item(item=item, save=False)
            self.effectors.append(node)

        elif isinstance(item, PathwayItem):

            node = Pathway.from_item(item=item, save=False)
            self.pathways.append(node)

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

        else:
            return item

        return item
