from typing import List

from protrend.models.node import Node
from protrend.models.version import Version
from protrend.utils.db_connection import DBSettings
from protrend.utils.node_importer import NodeImporter


class NeoPipeline:

    def __init__(self,
                 user_name: str,
                 password: str,
                 ip: str,
                 port: str,
                 db_name: str,
                 dbms: str,
                 import_folder: str,
                 version: str):

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
        if self._database is None:
            self._database = DBSettings(user_name=self._user_name,
                                        password=self._password,
                                        ip=self._ip,
                                        port=self._port,
                                        db_name=self._db_name,
                                        dbms=self._dbms)

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
