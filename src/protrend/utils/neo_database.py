import subprocess

from neomodel import config
from neomodel import db, clear_neo4j_database, install_all_labels, install_labels


class NeoDatabase:

    def __init__(self,
                 user_name: str,
                 password: str,
                 ip: str,
                 port: str,
                 db_name: str = '',
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

    def import_csv_data(self, *args):

        if not args:
            return

        cmd_arguments = [fr'{self.dbms}\bin\neo4j-admin',
                         'import',
                         '--database',
                         f'{self.db_name}',
                         '--multiline-fields=true']

        def sort_by_node(arg):
            if '--nodes=' in arg:
                return 0
            elif '--relationships=' in arg:
                return 1

            return 2

        sorted_args = sorted(args, key=sort_by_node)

        cmd_arguments.extend(sorted_args)

        return subprocess.run(cmd_arguments,
                              check=True,
                              text=True,
                              shell=True)
