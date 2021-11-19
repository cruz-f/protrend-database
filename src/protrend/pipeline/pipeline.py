from collections import defaultdict
from typing import List, Dict, Tuple, Type, TYPE_CHECKING

from tqdm import tqdm

from protrend.log import ProtrendLogger

if TYPE_CHECKING:
    from protrend.extract import Extractor
    from protrend.transform import Transformer
    from protrend.transform import Connector
    from protrend.load import Loader


def sort_by_order(transformer: 'Transformer'):
    return transformer.order


class Pipeline:
    default_extractors: Dict[Tuple[str, str], List[Type['Extractor']]] = defaultdict(list)
    default_transformers: Dict[Tuple[str, str], List[Type['Transformer']]] = defaultdict(list)
    default_connectors: Dict[Tuple[str, str], List[Type['Connector']]] = defaultdict(list)
    default_loaders: Dict[Tuple[str, str], List[Type['Loader']]] = defaultdict(list)

    @classmethod
    def register_extractor(cls, extractor: Type['Extractor'], source: str, version: str):
        cls.default_extractors[(source, version)].append(extractor)

    @classmethod
    def register_transformer(cls, transformer: Type['Transformer'], source: str, version: str):
        cls.default_transformers[(source, version)].append(transformer)

    @classmethod
    def register_connector(cls, connector: Type['Connector'], source: str, version: str):
        cls.default_connectors[(source, version)].append(connector)

    @classmethod
    def register_loader(cls, loader: Type['Loader'], source: str, version: str):
        cls.default_loaders[(source, version)].append(loader)

    def __init__(self,
                 extractors: List['Extractor'] = None,
                 transformers: List['Transformer'] = None,
                 connectors: List['Connector'] = None,
                 loaders: List['Loader'] = None):

        """
        Pipeline is responsible for managing transformers, connectors and loaders.
        It sends instructions to the transformers on how to transform a file(s) into a node or relationship
        Reading, processing, integrating and writing are executed by each transformer

        The director is also responsible for sending instructions to the loaders on how to load the transformed nodes
        into the database.

        :type extractors: List['Extractor']
        :param extractors: extractors to scrap data sources

        :type transformers: List['Transformer']
        :param transformers: transformers for processing staging area files into nodes

        :type connectors: List['Connector']
        :param connectors: connectors for processing staging area files into relationships

        :type loaders: List['Loader']
        :param loaders: loaders for converting data lake files into nodes and relationships
        """

        if not extractors:
            extractors = []

        if not transformers:
            transformers = []

        if not connectors:
            connectors = []

        if not loaders:
            loaders = []

        self._extractors = list(extractors)
        self._transformers = list(transformers)
        self._connectors = list(connectors)
        self._loaders = list(loaders)

    @property
    def extractors(self) -> List['Extractor']:
        return self._extractors

    @property
    def transformers(self) -> List['Transformer']:
        return sorted(self._transformers, key=sort_by_order, reverse=True)

    @property
    def connectors(self) -> List['Connector']:
        return self._connectors

    @property
    def loaders(self) -> List['Loader']:
        return self._loaders

    def extract(self):
        """
        Extracting data sources.

        :return:
        """

        extractors_info = ' '.join([extractor.source for extractor in self.extractors])
        ProtrendLogger.log.info(f'Pipeline call for extract using the spiders: {extractors_info}')

        for extractor in self.extractors:

            try:

                ProtrendLogger.log.info(f'Starting extractor: {extractor}')

                extractor.extract()

                ProtrendLogger.log.info(f'Finishing transformer: {extractor.source}')

            except:
                ProtrendLogger.log.exception(f'Exception occurred in extractor {extractor}')

    def transform(self):
        """
        Reading, processing, transforming, cleaning, integrating and writing one or more file types
        into multiple nodes.
        Transformation is performed step-wise according to the transformers order.

        :return:
        """

        transformers_info = ' '.join([transformer.node.node_name() for transformer in self.transformers])
        ProtrendLogger.log.info(f'Pipeline call for transform in the following transformers: {transformers_info}')

        for transformer in self.transformers:

            try:

                ProtrendLogger.log.info(f'Starting transformer: {transformer}')

                df = transformer.transform()
                transformer.integrate(df)
                transformer.write()

                ProtrendLogger.log.info(f'Finishing transformer: {transformer.node.node_name()} with stats: {df.shape}')

            except:
                ProtrendLogger.log.exception(f'Exception occurred in transformer {transformer}')

    def connect(self):
        """
        Reading, processing, transforming, cleaning, integrating and writing one or more file types
        into multiple relationships.

        :return:
        """

        connectors_info = ' '.join([f'{connector.from_node.node_name()}_{connector.to_node.node_name()}\n'
                                    for connector in self.connectors])
        ProtrendLogger.log.info(f'Pipeline call for connect in the following connectors: \n {connectors_info}')

        for connector in self.connectors:

            try:

                ProtrendLogger.log.info(f'Starting connector: {connector}')

                connector.connect()
                connector.write()

                ProtrendLogger.log.info(f'Finishing connector: '
                                        f'{connector.from_node.node_name()}_{connector.to_node.node_name()}')

            except:
                ProtrendLogger.log.exception(f'Exception occurred in connector {connector}')

    def load(self):
        """
        Loading multiple nodes and relationships. The nodes and relationships must be encoded into a csv file created
        by a transformer.
        Loading is performed step-wise according to the loaders order.

        :return:
        """

        loaders = ' '.join([f'{loader.__class__.__name__}' for loader in self.loaders])
        ProtrendLogger.log.info(f'Pipeline call for load in the following loaders: {loaders}')

        for loader in self.loaders:
            try:
                ProtrendLogger.log.info(f'Starting loader: {loader.__class__.__name__} with the following files:')
                for file_path in loader.load_stack:
                    ProtrendLogger.log.info(f'{file_path}')

                loader.load()

                ProtrendLogger.log.info(f'Finishing loader: {loader.__class__.__name__}')

            except:
                ProtrendLogger.log.exception(f'Exception occurred in loader {loader.__class__.__name__}')

    @classmethod
    def for_extraction(cls, source: str, version: str) -> 'Pipeline':
        extractors_types = cls.default_extractors.get((source, version), [])
        extractors = [extractor_type() for extractor_type in extractors_types]

        return cls(extractors=extractors)

    @classmethod
    def for_transformation(cls, source: str, version: str) -> 'Pipeline':
        transformers_types = cls.default_transformers.get((source, version), [])
        transformers = [transformer_type() for transformer_type in transformers_types]

        return cls(transformers=transformers)

    @classmethod
    def for_connection(cls, source: str, version: str) -> 'Pipeline':
        connectors_types = cls.default_connectors.get((source, version), [])
        connectors = [connector_type() for connector_type in connectors_types]

        return cls(connectors=connectors)

    @classmethod
    def for_loading(cls, source: str, version: str) -> 'Pipeline':
        loaders_types = cls.default_loaders.get((source, version), [])
        loaders = [loader_type() for loader_type in loaders_types]

        return cls(loaders=loaders)

    @classmethod
    def _for_tc(cls, source: str, version: str):
        transformers_types = cls.default_transformers.get((source, version), [])
        transformers = [transformer_type() for transformer_type in transformers_types]

        connectors_types = cls.default_connectors.get((source, version), [])
        connectors = [connector_type() for connector_type in connectors_types]

        return transformers, connectors

    @classmethod
    def _for_tcl(cls, source: str, version: str):
        transformers, connectors = cls._for_tc(source=source, version=version)

        loaders_types = cls.default_loaders.get((source, version), [])
        loaders = [loader_type() for loader_type in loaders_types]

        return transformers, connectors, loaders

    @classmethod
    def for_tc(cls, source: str, version: str) -> 'Pipeline':
        transformers, connectors, loaders = cls._for_tcl(source=source, version=version)
        return cls(transformers=transformers, connectors=connectors)

    @classmethod
    def for_tcl(cls, source: str, version: str) -> 'Pipeline':
        transformers, connectors, loaders = cls._for_tcl(source=source, version=version)
        return cls(transformers=transformers, connectors=connectors, loaders=loaders)

    @classmethod
    def for_etcl(cls, source: str, version: str) -> 'Pipeline':
        extractors_types = cls.default_extractors.get((source, version), [])
        extractors = [extractor_type() for extractor_type in extractors_types]

        transformers, connectors, loaders = cls._for_tcl(source=source, version=version)

        return cls(extractors=extractors, transformers=transformers, connectors=connectors, loaders=loaders)
