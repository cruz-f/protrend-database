from protrend.load import Loader


class AbasyLoader(Loader, source='abasy', version='0.0.0', register=True):
    pass


class CollecTFLoader(Loader, source='collectf', version='0.0.1', register=True):
    pass


class CoryneRegNetLoader(Loader, source='coryneregnet', version='0.0.0', register=True):
    pass


class DBTBSLoader(Loader, source='dbtbs', version='0.0.4', register=True):
    pass


class LiteratureLoader(Loader, source='literature', version='0.0.0', register=True):
    pass


class OperonDBLoader(Loader, source='operondb', version='0.0.0', register=True):
    pass


class MotifLoader(Loader, source='motif', version='0.0.0', register=True):
    pass


class RegPreciseLoader(Loader, source='regprecise', version='0.0.0', register=True):
    pass


class RegulonDBLoader(Loader, source='regulondb', version='0.0.0', register=True):
    pass


class StandardizerLoader(Loader, source='standardizer', version='0.0.0', register=True):
    pass
