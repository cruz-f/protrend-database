import logging.config
from datetime import datetime

from protrend.utils import ROOT_PATH

CONFIG_FILE = ROOT_PATH.joinpath('log', 'log.conf')
CURRENT_TIME = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
LOG_PATH = ROOT_PATH.joinpath('log', f'protrendTL_{CURRENT_TIME}.log')
LOG_FILE = LOG_PATH.as_posix()


def _singleton(cls):
    instances = {}

    def get_instance():
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]

    return get_instance()


@_singleton
class ProtrendLogger:

    def __init__(self):
        self._config = CONFIG_FILE
        self._disable_existing_loggers = False
        self._log_file = LOG_FILE
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
# is a Type an not an instance
ProtrendLogger: ProtrendLogger
