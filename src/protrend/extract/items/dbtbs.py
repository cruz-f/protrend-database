from itemloaders.processors import TakeFirst, Join, MapCompose
from scrapy import Item, Field

from protrend.extract.processors.dbtbs import DBTBSProcessors


class DatabaseItem(Item):
    name = Field()
    url = Field()
    doi = Field()
    authors = Field()
    description = Field()
    version = Field()
    created = Field()


class TranscriptionFactorItem(Item):
    family = Field()
    domain = Field()
    domain_description = Field()

    name = Field(output_processor=TakeFirst())
    description = Field()
    url = Field(output_processor=TakeFirst())

    type = Field()
    swiss_prot = Field(input_processor=MapCompose(DBTBSProcessors.process_nd))
    subti_list = Field(input_processor=MapCompose(DBTBSProcessors.process_nd))
    consensus_sequence = Field(input_processor=MapCompose(DBTBSProcessors.process_nd))
    comment = Field(input_processor=MapCompose(DBTBSProcessors.process_nd))

    # relationships
    operon = Field()
    gene = Field()
    tfbs = Field()


class OperonItem(Item):
    name = Field()
    url = Field()

    evidence = Field(input_processor=MapCompose(DBTBSProcessors.process_nd))
    pubmed = Field(input_processor=MapCompose(DBTBSProcessors.process_pubmed, DBTBSProcessors.process_nd))
    comment = Field(input_processor=MapCompose(DBTBSProcessors.process_nd))

    # relationships
    tf = Field()
    gene = Field()
    tfbs = Field()


class GeneItem(Item):
    name = Field()
    url = Field()
    synonyms = Field()
    strand = Field()
    position = Field()
    function = Field()
    cog_id = Field()
    conversed_groups = Field()

    # relationships
    tf = Field()
    operon = Field(output_processor=TakeFirst())
    tfbs = Field()


class TFBSItem(Item):
    identifier = Field(output_processor=Join(separator='_'))
    url = Field()
    regulation = Field()
    location = Field(input_processor=MapCompose(DBTBSProcessors.process_nd))
    absolute_position = Field(input_processor=MapCompose(DBTBSProcessors.process_nd))
    sequence = Field(input_processor=MapCompose(DBTBSProcessors.process_nd), output_processor=Join(separator=''))
    pubmed = Field(input_processor=MapCompose(DBTBSProcessors.process_pubmed, DBTBSProcessors.process_nd))

    # relationships
    tf = Field()
    operon = Field()
    gene = Field()
