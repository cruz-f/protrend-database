# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface
from typing import List

from protrend.extract.databases import RegPreciseDB
from protrend.extract.items import (TaxonomyItem, GenomeItem, TranscriptionFactorItem, RegulogItem,
                                    TranscriptionFactorFamilyItem, RNAFamilyItem, EffectorItem, PathwayItem,
                                    RegulonItem, OperonItem, GeneItem, TFBSItem)
from protrend.extract.utils import NodeRelationshipMap
from protrend.models.regprecise import (RegPreciseVersion, Taxonomy, TranscriptionFactor, TranscriptionFactorFamily,
                                        RNAFamily, Effector, Pathway, Genome, Regulog, Regulon, Operon, Gene, TFBS)
from protrend.models.regprecise import Database as RegPreciseDatabase
from protrend.models.version import Version
from protrend.models.node import Node
from protrend.utils.db_connection import DBSettings


class NeoPipeline:

    def __init__(self,
                 user_name: str = None,
                 password: str = None,
                 ip: str = None,
                 port: str = None,
                 version: str = None,
                 clear_version: bool = False,
                 clear_db: bool = False,
                 clear_schema: bool = False):

        if not user_name:
            user_name = 'db'

        if not password:
            password = 'db'

        if not ip:
            ip = 'localhost'

        if not port:
            port = '7687'

        self._user_name = user_name
        self._password = password
        self._ip = ip
        self._port = port

        self._database = None
        self._database_node = None

        self._version = version

        self.clear_version = clear_version
        self.clear_db = clear_db
        self.clear_schema = clear_schema

        self._relationships = []

    @classmethod
    def from_crawler(cls, crawler):

        user_name = crawler.settings.get('user_name')
        password = crawler.settings.get('password')
        ip = crawler.settings.get('ip')
        port = crawler.settings.get('port')
        version = crawler.settings.get('version')
        clear_version = crawler.settings.get('clear_version', False)
        clear_db = crawler.settings.get('clear_db', False)
        clear_schema = crawler.settings.get('clear_schema', False)

        return cls(user_name=user_name,
                   password=password,
                   ip=ip,
                   port=port,
                   version=version,
                   clear_version=clear_version,
                   clear_db=clear_db,
                   clear_schema=clear_schema)

    @property
    def relationships(self) -> List[NodeRelationshipMap]:
        return self._relationships

    @property
    def version(self) -> Version:

        return self._version

    @property
    def database(self) -> DBSettings:

        return self._database

    @property
    def database_node(self) -> Node:

        return self._database_node

    def open_spider(self, spider):
        pass

    def close_spider(self, spider):
        pass

    def process_item(self, item, spider):
        return item


class RegPrecisePipeline(NeoPipeline):

    def __init__(self,
                 user_name: str = None,
                 password: str = None,
                 ip: str = None,
                 port: str = None,
                 version: str = None,
                 clear_version: bool = False,
                 clear_db: bool = False,
                 clear_schema: bool = False):

        if not user_name:
            user_name = 'regprecise'

        if not password:
            password = 'regprecise'

        if not ip:
            ip = 'localhost'

        if not port:
            port = '7687'

        super().__init__(user_name=user_name,
                         password=password,
                         ip=ip,
                         port=port,
                         version=version,
                         clear_version=clear_version,
                         clear_db=clear_db,
                         clear_schema=clear_schema)

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

        for relationship in self.relationships:
            relationship.connect(self.version)

    def process_item(self, item, spider):

        if isinstance(item, TaxonomyItem):

            node, relationships = self.process_taxonomy_item(item=item, spider=spider)

        elif isinstance(item, GenomeItem):

            node, relationships = self.process_genome_item(item=item, spider=spider)

        elif isinstance(item, TranscriptionFactorItem):

            node, relationships = self.process_tf_item(item=item, spider=spider)

        elif isinstance(item, RegulogItem):

            node, relationships = self.process_regulog_item(item=item, spider=spider)

        elif isinstance(item, TranscriptionFactorFamilyItem):

            node, relationships = self.process_tffam_item(item=item, spider=spider)

        elif isinstance(item, RNAFamilyItem):

            node, relationships = self.process_rfam_item(item=item, spider=spider)

        elif isinstance(item, EffectorItem):

            node, relationships = self.process_effector_item(item=item, spider=spider)

        elif isinstance(item, PathwayItem):

            node, relationships = self.process_pathway_item(item=item, spider=spider)

        elif isinstance(item, RegulonItem):

            node, relationships = self.process_regulon_item(item=item, spider=spider)

        elif isinstance(item, OperonItem):

            node, relationships = self.process_operon_item(item=item, spider=spider)

        elif isinstance(item, GeneItem):

            node, relationships = self.process_gene_item(item=item, spider=spider)

        elif isinstance(item, TFBSItem):

            node, relationships = self.process_tfbs_item(item=item, spider=spider)

        else:
            return item

        node.connect_to_version(self.version)
        self.relationships.extend(relationships)

        return item

    def process_taxonomy_item(self, item, spider):
        attr = 'collection_id'

        node = Taxonomy.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:

            node = Taxonomy.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        node.database.connect(self.database_node)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='genome',
                                      to_node_cls=Genome,
                                      to_node_attr='genome_id',
                                      to=item['genome'])

        return node, [rel_map]

    def process_genome_item(self, item, spider):
        attr = 'genome_id'

        node = Genome.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Genome.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        rels = []

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='taxonomy',
                                      to_node_cls=Taxonomy,
                                      to_node_attr='collection_id',
                                      to=[item['taxonomy']])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='regulon',
                                      to_node_cls=Regulon,
                                      to_node_attr='regulon_id',
                                      to=item['regulon'])

        rels.append(rel_map)

        return node, rels

    def process_tf_item(self, item, spider):
        attr = 'collection_id'

        node = TranscriptionFactor.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = TranscriptionFactor.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        node.database.connect(self.database_node)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='regulog',
                                      to_node_cls=Regulog,
                                      to_node_attr='regulog_id',
                                      to=item['regulog'])

        return node, [rel_map]

    def process_regulog_item(self, item, spider):
        attr = 'regulog_id'

        node = Regulog.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Regulog.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        rels = []

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='transcription_factor',
                                      to_node_cls=TranscriptionFactor,
                                      to_node_attr='collection_id',
                                      to=[item['transcription_factor']])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='regulon',
                                      to_node_cls=Regulon,
                                      to_node_attr='regulon_id',
                                      to=item['regulon'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='taxonomy',
                                      to_node_cls=Taxonomy,
                                      to_node_attr='collection_id',
                                      to=item['taxonomy'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='tf_family',
                                      to_node_cls=TranscriptionFactorFamily,
                                      to_node_attr='tffamily_id',
                                      to=item['tf_family'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='rna_family',
                                      to_node_cls=RNAFamily,
                                      to_node_attr='riboswitch_id',
                                      to=item['rna_family'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='effector',
                                      to_node_cls=Effector,
                                      to_node_attr='effector_id',
                                      to=item['effector'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='pathway',
                                      to_node_cls=Pathway,
                                      to_node_attr='pathway_id',
                                      to=item['pathway'])

        rels.append(rel_map)

        return node, rels

    def process_tffam_item(self, item, spider):
        attr = 'tffamily_id'

        node = TranscriptionFactorFamily.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = TranscriptionFactorFamily.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        node.database.connect(self.database_node)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='regulog',
                                      to_node_cls=Regulog,
                                      to_node_attr='regulog_id',
                                      to=item['regulog'])

        return node, [rel_map]

    def process_rfam_item(self, item, spider):
        attr = 'riboswitch_id'

        node = RNAFamily.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = RNAFamily.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        node.database.connect(self.database_node)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='regulog',
                                      to_node_cls=Regulog,
                                      to_node_attr='regulog_id',
                                      to=item['regulog'])

        return node, [rel_map]

    def process_effector_item(self, item, spider):
        attr = 'effector_id'

        node = Effector.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Effector.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        node.database.connect(self.database_node)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='regulog',
                                      to_node_cls=Regulog,
                                      to_node_attr='regulog_id',
                                      to=item['regulog'])

        return node, [rel_map]

    def process_pathway_item(self, item, spider):
        attr = 'pathway_id'

        node = Pathway.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Pathway.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        node.database.connect(self.database_node)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='regulog',
                                      to_node_cls=Regulog,
                                      to_node_attr='regulog_id',
                                      to=item['regulog'])

        return node, [rel_map]

    def process_regulon_item(self, item, spider):
        attr = 'regulon_id'

        node = Regulon.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Regulon.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        rels = []

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='genome',
                                      to_node_cls=Genome,
                                      to_node_attr='genome_id',
                                      to=[item['genome']])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='operon',
                                      to_node_cls=Operon,
                                      to_node_attr='operon_id',
                                      to=item['operon'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='gene',
                                      to_node_cls=Gene,
                                      to_node_attr='locus_tag',
                                      to=item['gene'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='tfbs',
                                      to_node_cls=TFBS,
                                      to_node_attr='tfbs_id',
                                      to=item['tfbs'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='regulog',
                                      to_node_cls=Regulog,
                                      to_node_attr='regulog_id',
                                      to=item['regulog'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='taxonomy',
                                      to_node_cls=Taxonomy,
                                      to_node_attr='collection_id',
                                      to=item['taxonomy'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='transcription_factor',
                                      to_node_cls=TranscriptionFactor,
                                      to_node_attr='collection_id',
                                      to=item['transcription_factor'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='tf_family',
                                      to_node_cls=TranscriptionFactorFamily,
                                      to_node_attr='tffamily_id',
                                      to=item['tf_family'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='rna_family',
                                      to_node_cls=RNAFamily,
                                      to_node_attr='riboswitch_id',
                                      to=item['rna_family'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='effector',
                                      to_node_cls=Effector,
                                      to_node_attr='effector_id',
                                      to=item['effector'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='pathway',
                                      to_node_cls=Pathway,
                                      to_node_attr='pathway_id',
                                      to=item['pathway'])

        rels.append(rel_map)

        return node, rels

    def process_operon_item(self, item, spider):
        attr = 'operon_id'

        node = Operon.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Operon.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        rels = []

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='regulon',
                                      to_node_cls=Regulon,
                                      to_node_attr='regulon_id',
                                      to=[item['regulon']])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='gene',
                                      to_node_cls=Gene,
                                      to_node_attr='locus_tag',
                                      to=item['gene'])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='tfbs',
                                      to_node_cls=TFBS,
                                      to_node_attr='tfbs_id',
                                      to=item['tfbs'])

        rels.append(rel_map)

        return node, rels

    def process_gene_item(self, item, spider):
        attr = 'locus_tag'

        node = Gene.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Gene.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        rels = []

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='regulon',
                                      to_node_cls=Regulon,
                                      to_node_attr='regulon_id',
                                      to=[item['regulon']])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='operon',
                                      to_node_cls=Operon,
                                      to_node_attr='operon_id',
                                      to=[item['operon']])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='tfbs',
                                      to_node_cls=TFBS,
                                      to_node_attr='tfbs_id',
                                      to=item['tfbs'])

        rels.append(rel_map)

        return node, rels

    def process_tfbs_item(self, item, spider):
        attr = 'tfbs_id'

        node = TFBS.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = TFBS.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        rels = []

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='regulon',
                                      to_node_cls=Regulon,
                                      to_node_attr='regulon_id',
                                      to=[item['regulon']])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='operon',
                                      to_node_cls=Operon,
                                      to_node_attr='operon_id',
                                      to=[item['operon']])

        rels.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      relationship_name='gene',
                                      to_node_cls=TFBS,
                                      to_node_attr='locus_tag',
                                      to=item['gene'])

        rels.append(rel_map)

        return node, rels
