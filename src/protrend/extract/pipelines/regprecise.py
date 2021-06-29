from protrend.extract.databases import RegPreciseDB
from protrend.extract.items.regprecise import (TaxonomyItem, GenomeItem, TranscriptionFactorItem, RegulogItem,
                                               TranscriptionFactorFamilyItem, RNAFamilyItem, EffectorItem, PathwayItem,
                                               RegulonItem, OperonItem, GeneItem, TFBSItem)
from protrend.extract.pipelines.neo_pipeline import NeoPipeline
from protrend.models.regprecise import Database as RegPreciseDatabase
from protrend.models.regprecise import (RegPreciseVersion, Taxonomy, TranscriptionFactor, TranscriptionFactorFamily,
                                        RNAFamily, Effector, Pathway, Genome, Regulog, Regulon, Operon, Gene, TFBS)
from protrend.utils.node_importer import NodeImporter


class RegPrecisePipeline(NeoPipeline):

    def __init__(self,
                 user_name: str = None,
                 password: str = None,
                 *args,
                 **kwargs):

        if not user_name:
            user_name = 'regprecise'

        if not password:
            password = 'regprecise'

        super().__init__(user_name=user_name,
                         password=password,
                         *args,
                         **kwargs)

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
    def database(self) -> RegPreciseDB:

        if self._database is None:
            self._database = RegPreciseDB(user_name=self._user_name,
                                          password=self._password,
                                          ip=self._ip,
                                          port=self._port,
                                          db_name=self._db_name,
                                          dbms=self._dbms)

        return self._database

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
    def database_node(self) -> RegPreciseDatabase:

        if self._database_node is None:

            node = RegPreciseDatabase.get_by_version(attr='name', value='regprecise', version=self.version)

            if node is None:
                node = RegPreciseDatabase().save()

            self._database_node = node

        return self._database_node

    def open_spider(self, spider):

        self.database.connect()

        if self.clear_schema:
            self.database.clear_db(clear_constraints=True, clear_indexes=True)
            self.database.install_all_labels()

        elif self.clear_db:
            self.database.clear_db()

        elif self.clear_version:

            children = self.version.get_children()

            for node in children:
                node.delete()

        self.database_node.connect_to_version(self.version)

    def close_spider(self, spider):

        arguments = []

        for importer in self.importers:
            importer.build_imports()
            arguments.append(importer.args)

        self.database.import_csv_data(arguments=arguments)

    def process_item(self, item, spider):

        if isinstance(item, TaxonomyItem):

            node = Taxonomy.from_item(item=item, version=self.version.name, save=False)
            self.taxa.append(node)

        elif isinstance(item, GenomeItem):

            node = Genome.from_item(item=item, version=self.version.name, save=False)
            self.genomes.append(node)

        elif isinstance(item, TranscriptionFactorItem):

            node = TranscriptionFactor.from_item(item=item, version=self.version.name, save=False)
            self.tfs.append(node)

        elif isinstance(item, RegulogItem):

            node = Regulog.from_item(item=item, version=self.version.name, save=False)
            self.regulogs.append(node)

        elif isinstance(item, TranscriptionFactorFamilyItem):

            node = TranscriptionFactorFamily.from_item(item=item, version=self.version.name, save=False)
            self.tffams.append(node)

        elif isinstance(item, RNAFamilyItem):

            node = RNAFamily.from_item(item=item, version=self.version.name, save=False)
            self.rnafams.append(node)

        elif isinstance(item, EffectorItem):

            node = Effector.from_item(item=item, version=self.version.name, save=False)
            self.effectors.append(node)

        elif isinstance(item, PathwayItem):

            node = Pathway.from_item(item=item, version=self.version.name, save=False)
            self.pathways.append(node)

        elif isinstance(item, RegulonItem):

            node = Regulon.from_item(item=item, version=self.version.name, save=False)
            self.regulons.append(node)

        elif isinstance(item, OperonItem):

            node = Operon.from_item(item=item, version=self.version.name, save=False)
            self.operons.append(node)

        elif isinstance(item, GeneItem):

            node = Gene.from_item(item=item, version=self.version.name, save=False)
            self.genes.append(node)

        elif isinstance(item, TFBSItem):

            node = TFBS.from_item(item=item, version=self.version.name, save=False)
            self.tfbs.append(node)

        else:
            return item

        return item
