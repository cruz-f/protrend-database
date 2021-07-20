class BioAPI:

    def __init__(self):
        self._record = {}

    @property
    def record(self):
        return self._record

    @record.setter
    def record(self, value):
        if value:
            self._record = value

    def is_empty(self):
        if self.record:
            return True
        return False

    def fetch(self):
        pass
