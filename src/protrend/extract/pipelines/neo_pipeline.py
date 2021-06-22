from typing import List

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