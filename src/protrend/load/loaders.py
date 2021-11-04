from protrend.load import Loader


class AbasyLoader(Loader, source='abasy', version='0.0.0', register=True):
    pass


class CollectfLoader(Loader, source='collectf', version='0.0.1', register=True):
    pass


class CoryneRegNetLoader(Loader, source='coryneregnet', version='0.0.0', register=True):
    pass


class DBTBSLoader(Loader, source='dbtbs', version='0.0.3', register=True):
    pass


class LiteratureLoader(Loader, source='literature', version='0.0.0', register=True):
    pass


class RegPreciseLoader(Loader, source='regprecise', version='0.0.0', register=True):
    pass


class RegulondbLoader(Loader, source='regulondb', version='0.0.0', register=True):
    pass
