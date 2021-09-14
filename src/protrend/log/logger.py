import logging.config
import os
from datetime import datetime

from protrend.utils import ROOT_PATH

CONFIG_FILE = ROOT_PATH.joinpath('log', 'log.conf')
CURRENT_TIME = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
LOG_FILE = ROOT_PATH.joinpath('log', f'protrendTL_{CURRENT_TIME}.log')
LOG_FILE = os.fspath(LOG_FILE)


def singleton(cls):
    instances = {}

    def get_instance():
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]

    return get_instance()


@singleton
class Logger:

    def __init__(self):
        logging.config.fileConfig(fname=CONFIG_FILE, disable_existing_loggers=False,
                                  defaults={'logfilename': LOG_FILE})
        self.log = logging.getLogger('protrendTL')
        self.log.info('Started logger')
