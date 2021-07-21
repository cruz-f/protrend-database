class BioAPI:

    def __init__(self, identifier: str):

        if not identifier:
            identifier = ''

        self._identifier = identifier
        self._record = {}

    @property
    def identifier(self) -> str:
        return self._identifier

    @property
    def record(self) -> dict:
        return self._record

    @record.setter
    def record(self, value):
        if value:
            self._record = value

    def is_empty(self) -> bool:
        if self.record:
            return True
        return False

    def fetch(self, *args, **kwargs):
        pass
