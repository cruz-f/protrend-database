import os
from pathlib import Path
from typing import List

from protrend.models.node import Node
from protrend.models.version import Version
from protrend.utils.db_connection import DBSettings
from protrend.utils.node_importer import NodeImporter


class NeoPipeline:

    def __init__(self,
                 user_name: str = None,
                 password: str = None,
                 ip: str = None,
                 port: str = None,
                 db_name: str = None,
                 dbms: str = None,
                 import_folder: str = None,
                 version: str = None):

        if not user_name:
            user_name = 'db'

        if not password:
            password = 'db'

        if not ip:
            ip = 'localhost'

        if not port:
            port = '7687'

        if not db_name:
            db_name = 'neo4j'

        if not dbms:
            dbms = ''

        if not import_folder:
            import_path = Path(os.path.abspath(__file__))
            import_folder = os.path.join(str(import_path.parent.absolute()), 'import')

        self._user_name = user_name
        self._password = password
        self._ip = ip
        self._port = port
        self._db_name = db_name
        self._dbms = dbms

        self._import_folder = import_folder

        self._database = None
        self._database_node = None

        self._version = version

    @classmethod
    def from_crawler(cls, crawler):

        user_name = crawler.settings.get('user_name')
        password = crawler.settings.get('password')
        ip = crawler.settings.get('ip')
        port = crawler.settings.get('port')
        db_name = crawler.settings.get('db_name')
        dbms = crawler.settings.get('dbms')
        import_folder = crawler.settings.get('import_folder')
        version = crawler.settings.get('version')

        return cls(user_name=user_name,
                   password=password,
                   ip=ip,
                   port=port,
                   db_name=db_name,
                   dbms=dbms,
                   import_folder=import_folder,
                   version=version)

    @property
    def importers(self) -> List[NodeImporter]:
        return []

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
