from protrend.extract import Extractor


class RegPreciseExtractor(Extractor, source='regprecise', version='0.0.0', register=True):
    pass


class CollectfExtractor(Extractor, source='collectf', version='0.0.1', register=True):
    pass


class DBTBSExtractor(Extractor, source='dbtbs', version='0.0.3', register=True):
    pass
