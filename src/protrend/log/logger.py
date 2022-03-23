import logging.config
from datetime import datetime

from protrend.utils import Settings, singleton


@singleton
class ProtrendLogger:

    def __init__(self):
        self._config = Settings.log_conf_file
        self._disable_existing_loggers = False

        now = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
        self._log_file = Settings.log_working_directory.joinpath(f'protrendTL_{now}.log').as_posix()

        self._log = None
        self._started = False

    @property
    def config(self):
        return self._config

    @config.setter
    def config(self, value):

        if self._started:
            raise ValueError('ProtrendLogger configuration cannot be changed after starting the logger')

        self._config = value

    @property
    def disable_existing_loggers(self):
        return self._disable_existing_loggers

    @disable_existing_loggers.setter
    def disable_existing_loggers(self, value):

        if self._started:
            raise ValueError('ProtrendLogger configuration cannot be changed after starting the logger')

        self._disable_existing_loggers = value

    @property
    def log_file(self):
        return self._log_file

    @log_file.setter
    def log_file(self, value):

        if self._started:
            raise ValueError('ProtrendLogger configuration cannot be changed after starting the logger')

        self._log_file = value

    @property
    def defaults(self):
        return {'logfilename': self.log_file}

    @property
    def log(self) -> logging.Logger:
        return self._log

    def start_logger(self):
        logging.config.fileConfig(fname=self.config,
                                  disable_existing_loggers=self.disable_existing_loggers,
                                  defaults=self.defaults)
        self._log = logging.getLogger('protrendTL')
        self.log.info('Started logger')
        self._started = True


# Due to the singleton pattern, code inspections, lints and stubs will raise warnings that ProtrendLogger
# is a Type and not an instance
ProtrendLogger: ProtrendLogger
