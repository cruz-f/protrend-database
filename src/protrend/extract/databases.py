from protrend.utils.db_connection import DBSettings


class RegPreciseDB(DBSettings):

    def __init__(self,
                 user_name: str = 'regprecise',
                 password: str = 'regprecise',
                 ip: str = 'localhost',
                 port: str = '7687'):

        super().__init__(user_name, password, ip, port)