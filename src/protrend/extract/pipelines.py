# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface
from protrend.extract.databases import RegPreciseDB
from protrend.models.regprecise import Version as RegPreciseVersion
from protrend.models.version import VersionNode
from protrend.utils.db_connection import DBSettings


class NeoPipeline:

    def __init__(self,
                 user_name: str = None,
                 password: str = None,
                 ip: str = None,
                 port: str = None,
                 version: str = None,
                 clear_version: bool = False,
                 clear_sa: bool = False,
                 clear_sa_schema: bool = False):

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

        self._version = version

        self.clear_version = clear_version
        self.clear_sa = clear_sa
        self.clear_sa_schema = clear_sa_schema

    @classmethod
    def from_crawler(cls, crawler):

        user_name = crawler.settings.get('user_name')
        password = crawler.settings.get('password')
        ip = crawler.settings.get('ip')
        port = crawler.settings.get('port')
        version = crawler.settings.get('version')
        clear_version = crawler.settings.get('clear_version', False)
        clear_sa = crawler.settings.get('clear_sa', False)
        clear_sa_schema = crawler.settings.get('clear_sa_schema', False)

        return cls(user_name=user_name,
                   password=password,
                   ip=ip,
                   port=port,
                   version=version,
                   clear_version=clear_version,
                   clear_sa=clear_sa,
                   clear_sa_schema=clear_sa_schema)

    @property
    def version(self) -> VersionNode:

        return self._version

    @property
    def database(self) -> DBSettings:

        return self._database

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
                 clear_sa: bool = False,
                 clear_sa_schema: bool = False):

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
                         clear_sa=clear_sa,
                         clear_sa_schema=clear_sa_schema)

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
    def database(self) -> RegPreciseDB:

        if self._database is None:
            self._database = RegPreciseDB(user_name=self._user_name,
                                          password=self._password,
                                          ip=self._ip,
                                          port=self._port)

        return self._database

    def open_spider(self, spider):

        self.database.connect()

        if self.clear_sa_schema:
            self.database.clear_db(clear_constraints=True, clear_indexes=True)
            self.database.install_all_labels()

        elif self.clear_sa:
            self.database.clear_db()

        elif self.clear_version:

            children = self.version.get_children()

            for node in children:
                node.delete()

    def process_item(self, item, spider):

        return item
