from protrend.utils.db_connection import DBSettings


class RegPreciseDB(DBSettings):

    def __init__(self,
                 user_name: str = 'regprecise',
                 password: str = 'regprecise',
                 *args,
                 **kwargs):

        super().__init__(user_name, password, *args, **kwargs)
