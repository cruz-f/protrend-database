import json
import os
from datetime import datetime

from protrend.utils import singleton, Settings


@singleton
class ProtrendReporter:

    def __init__(self):
        now = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
        self._report_file = Settings.report_working_directory.joinpath(f'protrend_report_{now}.json')

        self._report = {}

    @property
    def report_file(self):
        return self._report_file

    @property
    def report(self):
        return self._report

    def report_objects(self,
                       source: str,
                       version: str,
                       system: str,
                       label: str,
                       objects: int,
                       properties: int):
        self._report[(source, version, system, label)] = {
            'source': source,
            'version': version,
            'system': system,
            'label': label,
            'objects': objects,
            'properties': properties}
        return

    def report_relationships(self,
                             source: str,
                             version: str,
                             system: str,
                             source_label: str,
                             target_label: str,
                             relationships: int):
        self._report[(source, version, system, source_label, target_label)] = {
            'source': source,
            'version': version,
            'system': system,
            'source_label': source_label,
            'target_label': target_label,
            'relationships': relationships}
        return

    def write_report(self):

        if not Settings.report_working_directory.exists():
            os.makedirs(Settings.report_working_directory)

        report = {'_'.join(key): value for key, value in self.report.items()}

        with open(self.report_file, 'w') as f:
            json.dump(report, f)

        # once the report is written, clear it
        self._report = {}
        return True


# Due to the singleton pattern, code inspections, lints and stubs will raise warnings that ProtrendLogger
# is a Type and not an instance
ProtrendReporter: ProtrendReporter
