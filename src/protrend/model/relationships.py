from neomodel import StructuredRel, DateTimeProperty, StringProperty

BASE_REL_TYPE = 'HAS'
SOURCE_REL_TYPE = 'OWNER'


class BaseRelationship(StructuredRel):
    # base
    created = DateTimeProperty(default_now=True)
    updated = DateTimeProperty(default_now=True)


class SourceRelationship(StructuredRel):
    # base
    created = DateTimeProperty(default_now=True)
    updated = DateTimeProperty(default_now=True)

    # properties
    key = StringProperty()
    url = StringProperty()
    external_identifier = StringProperty()
