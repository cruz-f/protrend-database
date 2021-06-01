from neomodel import db, clear_neo4j_database, install_all_labels, install_labels
from neomodel import config


class DBSettings:

    def __init__(self,
                 user_name: str = 'neo4j',
                 password: str = 'neo4j',
                 ip: str = 'localhost',
                 port: str = '7687'):

        self.user_name: str = user_name
        self.password: str = password
        self.ip: str = ip
        self.port: str = port

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
