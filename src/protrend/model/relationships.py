from neomodel import StructuredRel, DateTimeProperty, StringProperty

REL_TYPE = 'HAS'


class SourceRelationship(StructuredRel):
    # base
    created = DateTimeProperty(default_now=True)
    updated = DateTimeProperty(default_now=True)

    # properties
    key = StringProperty()
    url = StringProperty()
    external_identifier = StringProperty()
