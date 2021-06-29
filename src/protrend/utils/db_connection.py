import subprocess
from collections import namedtuple
from typing import Union, List, Set, Tuple

from neomodel import db, clear_neo4j_database, install_all_labels, install_labels
from neomodel import config

NeoImportEntity = namedtuple('NeoImportEntity', field_names=('label', 'csv_file'))

NeoImportsTyping = Union[List[NeoImportEntity], Set[NeoImportEntity], Tuple[NeoImportEntity]]


class DBSettings:

    def __init__(self,
                 user_name: str = 'neo4j',
                 password: str = 'neo4j',
                 ip: str = 'localhost',
                 port: str = '7687',
                 db_name: str = 'neo4j',
                 dbms: str = ''):

        self.user_name: str = user_name
        self.password: str = password
        self.ip: str = ip
        self.port: str = port
        self.db_name: str = db_name
        self.dbms: str = dbms

    @property
    def db(self) -> db:
        return db

    @property
    def bolt_url(self) -> str:
        return f'bolt://{self.user_name}:{self.password}@{self.ip}:{self.port}'

    def connect(self):
        db.set_connection(self.bolt_url)

    @staticmethod
    def auto_install_labels():
        config.AUTO_INSTALL_LABELS = True

    @staticmethod
    def install_labels(cls):
        install_labels(cls)

    @staticmethod
    def install_all_labels():
        install_all_labels()

    @staticmethod
    def clear_db(clear_constraints=False, clear_indexes=False):
        clear_neo4j_database(db, clear_constraints, clear_indexes)

    def import_csv_data(self, arguments):

        if not arguments:
            return

        cmd_arguments = [fr"{self.dbms}\bin\neo4j-admin",
                         "import",
                         "--database",
                         f"{self.db_name}"]

        cmd_arguments.extend(arguments)

        return subprocess.run(cmd_arguments,
                              check=True,
                              text=True,
                              shell=True)
