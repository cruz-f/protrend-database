import pandas as pd

from protrend.extract.databases import RegPreciseDB
from protrend.extract.items.regprecise import (TaxonomyItem, GenomeItem, TranscriptionFactorItem, RegulogItem,
                                               TranscriptionFactorFamilyItem, RNAFamilyItem, EffectorItem, PathwayItem,
                                               RegulonItem, OperonItem, GeneItem, TFBSItem)
from protrend.extract.pipelines.neo_pipeline import NeoPipeline
from protrend.models.regprecise import Database as RegPreciseDatabase
from protrend.models.regprecise import (RegPreciseVersion, Taxonomy, TranscriptionFactor, TranscriptionFactorFamily,
                                        RNAFamily, Effector, Pathway, Genome, Regulog, Regulon, Operon, Gene, TFBS)


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

        self.taxa = []
        self.genomes = []
        self.tfs = []
        self.regulogs = []
        self.tffams = []
        self.rnafams = []
        self.effectors = []
        self.pathways = []
        self.regulons = []
        self.operons = []
        self.genes = []
        self.tfbs = []

    @property
    def database(self) -> RegPreciseDB:

        if self._database is None:
            self._database = RegPreciseDB(user_name=self._user_name,
                                          password=self._password,
                                          ip=self._ip,
                                          port=self._port)

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

        nodes_containers = [(self.taxa, Taxonomy),
                            (self.genomes, Genome),
                            (self.tfs, TranscriptionFactor),
                            (self.regulogs, Regulog),
                            (self.tffams, TranscriptionFactorFamily),
                            (self.rnafams, RNAFamily),
                            (self.effectors, Effector),
                            (self.pathways, Pathway),
                            (self.regulons, Regulon),
                            (self.operons, Operon),
                            (self.genes, Gene),
                            (self.tfbs, TFBS)]

    def process_nodes_to_df(self, nodes, node_cls):
        series = [node.to_series() for node in nodes]
        df = pd.DataFrame(series, columns=node_cls.cls_keys())
        df.rename(columns={node_cls.property_as_id: f'{node_cls.property_as_id}:ID({node_cls.cls_name()})'},
                  inplace=True)
        return df

    def process_relationships_to_df(self, nodes, node_cls):
        pass

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
