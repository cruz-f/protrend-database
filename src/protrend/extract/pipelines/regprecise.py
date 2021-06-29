from protrend.extract.databases import RegPreciseDB
from protrend.extract.items.regprecise import (TaxonomyItem, GenomeItem, TranscriptionFactorItem, RegulogItem,
                                               TranscriptionFactorFamilyItem, RNAFamilyItem, EffectorItem, PathwayItem,
                                               RegulonItem, OperonItem, GeneItem, TFBSItem)
from protrend.extract.pipelines.neo_pipeline import NeoPipeline
from protrend.extract.utils import NodeRelationshipMap
from protrend.models.regprecise import (RegPreciseVersion, Taxonomy, TranscriptionFactor, TranscriptionFactorFamily,
                                        RNAFamily, Effector, Pathway, Genome, Regulog, Regulon, Operon, Gene, TFBS)
from protrend.models.regprecise import Database as RegPreciseDatabase


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

            node = Taxonomy.from_item(item=item)

        else:
            node = node.update(item=item)

        node.database.connect(self.database_node)

        relationships = [NodeRelationshipMap(node=node,
                                             item=item,
                                             relationship_name='genome',
                                             to_node_cls=Genome,
                                             to_node_attr='genome_id')]

        return node, relationships

    def process_genome_item(self, item, spider):
        attr = 'genome_id'

        node = Genome.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Genome.from_item(item=item)

        else:
            node = node.update(item=item)

        relationships = []

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='taxonomy',
                                      to_node_cls=Taxonomy,
                                      to_node_attr='collection_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='regulon',
                                      to_node_cls=Regulon,
                                      to_node_attr='regulon_id')

        relationships.append(rel_map)

        return node, relationships

    def process_tf_item(self, item, spider):
        attr = 'collection_id'

        node = TranscriptionFactor.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = TranscriptionFactor.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        node.database.connect(self.database_node)

        relationships = [NodeRelationshipMap(node=node,
                                             item=item,
                                             relationship_name='regulog',
                                             to_node_cls=Regulog,
                                             to_node_attr='regulog_id')]

        return node, relationships

    def process_regulog_item(self, item, spider):
        attr = 'regulog_id'

        node = Regulog.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Regulog.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        relationships = []

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='transcription_factor',
                                      to_node_cls=TranscriptionFactor,
                                      to_node_attr='collection_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='rna_family',
                                      to_node_cls=RNAFamily,
                                      to_node_attr='riboswitch_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='regulon',
                                      to_node_cls=Regulon,
                                      to_node_attr='regulon_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='taxonomy',
                                      to_node_cls=Taxonomy,
                                      to_node_attr='collection_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='tf_family',
                                      to_node_cls=TranscriptionFactorFamily,
                                      to_node_attr='tffamily_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='effector',
                                      to_node_cls=Effector,
                                      to_node_attr='effector_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='pathway',
                                      to_node_cls=Pathway,
                                      to_node_attr='pathway_id')

        relationships.append(rel_map)

        return node, relationships

    def process_tffam_item(self, item, spider):
        attr = 'tffamily_id'

        node = TranscriptionFactorFamily.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = TranscriptionFactorFamily.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        node.database.connect(self.database_node)

        relationships = [NodeRelationshipMap(node=node,
                                             item=item,
                                             relationship_name='regulog',
                                             to_node_cls=Regulog,
                                             to_node_attr='regulog_id')]

        return node, relationships

    def process_rfam_item(self, item, spider):
        attr = 'riboswitch_id'

        node = RNAFamily.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = RNAFamily.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        node.database.connect(self.database_node)

        relationships = [NodeRelationshipMap(node=node,
                                             item=item,
                                             relationship_name='regulog',
                                             to_node_cls=Regulog,
                                             to_node_attr='regulog_id')]

        return node, relationships

    def process_effector_item(self, item, spider):
        attr = 'effector_id'

        node = Effector.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Effector.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        node.database.connect(self.database_node)

        relationships = [NodeRelationshipMap(node=node,
                                             item=item,
                                             relationship_name='regulog',
                                             to_node_cls=Regulog,
                                             to_node_attr='regulog_id')]

        return node, relationships

    def process_pathway_item(self, item, spider):
        attr = 'pathway_id'

        node = Pathway.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Pathway.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        node.database.connect(self.database_node)

        relationships = [NodeRelationshipMap(node=node,
                                             item=item,
                                             relationship_name='regulog',
                                             to_node_cls=Regulog,
                                             to_node_attr='regulog_id')]

        return node, relationships

    def process_regulon_item(self, item, spider):
        attr = 'regulon_id'

        node = Regulon.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Regulon.from_item(item=item)

        else:
            node = node.update(item=item)

        relationships = []

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='genome',
                                      to_node_cls=Genome,
                                      to_node_attr='genome_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='operon',
                                      to_node_cls=Operon,
                                      to_node_attr='operon_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='gene',
                                      to_node_cls=Gene,
                                      to_node_attr='locus_tag')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='tfbs',
                                      to_node_cls=TFBS,
                                      to_node_attr='tfbs_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='regulog',
                                      to_node_cls=Regulog,
                                      to_node_attr='regulog_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='taxonomy',
                                      to_node_cls=Taxonomy,
                                      to_node_attr='collection_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='transcription_factor',
                                      to_node_cls=TranscriptionFactor,
                                      to_node_attr='collection_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='tf_family',
                                      to_node_cls=TranscriptionFactorFamily,
                                      to_node_attr='tffamily_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='rna_family',
                                      to_node_cls=RNAFamily,
                                      to_node_attr='riboswitch_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='effector',
                                      to_node_cls=Effector,
                                      to_node_attr='effector_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='pathway',
                                      to_node_cls=Pathway,
                                      to_node_attr='pathway_id')

        relationships.append(rel_map)

        return node, relationships

    def process_operon_item(self, item, spider):
        attr = 'operon_id'

        node = Operon.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Operon.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        relationships = []

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='regulon',
                                      to_node_cls=Regulon,
                                      to_node_attr='regulon_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='gene',
                                      to_node_cls=Gene,
                                      to_node_attr='locus_tag')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='tfbs',
                                      to_node_cls=TFBS,
                                      to_node_attr='tfbs_id')

        relationships.append(rel_map)

        return node, relationships

    def process_gene_item(self, item, spider):
        attr = 'locus_tag'

        node = Gene.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = Gene.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        relationships = []

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='regulon',
                                      to_node_cls=Regulon,
                                      to_node_attr='regulon_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='operon',
                                      to_node_cls=Operon,
                                      to_node_attr='operon_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='tfbs',
                                      to_node_cls=TFBS,
                                      to_node_attr='tfbs_id')

        relationships.append(rel_map)

        return node, relationships

    def process_tfbs_item(self, item, spider):
        attr = 'tfbs_id'

        node = TFBS.get_by_version(attr=attr, value=item[attr], version=self.version)

        if node is None:
            node = TFBS.from_item(item=item, save=True)

        else:
            node = node.update(item=item, save=True)

        relationships = []

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='regulon',
                                      to_node_cls=Regulon,
                                      to_node_attr='regulon_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='operon',
                                      to_node_cls=Operon,
                                      to_node_attr='operon_id')

        relationships.append(rel_map)

        rel_map = NodeRelationshipMap(node=node,
                                      item=item,
                                      relationship_name='gene',
                                      to_node_cls=Gene,
                                      to_node_attr='locus_tag')

        relationships.append(rel_map)

        return node, relationships
