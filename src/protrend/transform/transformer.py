import os
from abc import ABCMeta, abstractmethod
from functools import partial
from typing import List, Type, Callable, Dict

import pandas as pd

from protrend import GeneDTO, annotate_genes, OrganismDTO, annotate_organisms
from protrend.io import read_from_stack
from protrend.io.json import write_json_frame
from protrend.model.node import Node, protrend_id_decoder, protrend_id_encoder
from protrend.utils import SetList, Settings, WriteStack, DefaultProperty
from protrend.utils.processors import take_last, apply_processors, to_list_nan, regulatory_interaction_hash


def get_values(df, col):
    series = df.get(col)

    if series is None:
        return

    return series.to_list()


class AbstractTransformer(metaclass=ABCMeta):
    """
    Transformer interface.
    The following methods must be implemented to set up a transformer for each node
    """

    # --------------------------------------------------------
    # Transformer API
    # --------------------------------------------------------
    @abstractmethod
    def transform(self) -> pd.DataFrame:
        """
        The method responsible for transforming multiple pandas DataFrames into an annotated, cleaned and standardized
        pandas DataFrame ready to be integrated into structured nodes.

        Interface implementation with unknown signature.
        Concrete implementations are available at the transformer or transformer's sub-classes.

        :return:
        """

        pass

    @abstractmethod
    def integrate(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        The method responsible for integrating the transformed pandas DataFrame with the current state of the neo4j
        database using neomodel and Node class extension.

        Interface implementation with unknown signature.
        Concrete implementations are available at the transformer or transformer's sub-classes.

        :param df: transformed pandas DataFrame
        :return: it creates a new pandas DataFrame of the integrated data
        """
        pass

    def write(self):
        """
        The method responsible for writing the transformed, integrated and nodes pandas DataFrame
        to be loaded and connected in the neo4j database using neomodel and Node class extension.

        Interface implementation with unknown signature.
        Concrete implementations are available at the transformer or transformer's sub-classes.

        """
        pass


class Transformer(AbstractTransformer):
    """
    A transformer is responsible for reading, processing, transforming, integrating and writing node files.
    A transformer starts with data from the staging area and ends with structured nodes.

    This Transformer object has implemented several utilities for the hard-working de facto transformers to be created
    for each data source and node
    """
    source = DefaultProperty()
    version = DefaultProperty()
    node = DefaultProperty()
    order = DefaultProperty()

    default_transform_stack = {}
    columns = SetList()
    read_columns = SetList()

    def __init_subclass__(cls, **kwargs):

        source = kwargs.get('source')
        cls.source.set_default(cls, source)

        version = kwargs.get('version')
        cls.version.set_default(cls, version)

        node = kwargs.pop('node', None)
        cls.node.set_default(cls, node)

        order = kwargs.pop('order', None)
        cls.order.set_default(cls, order)

        register = kwargs.pop('register', False)

        if register:
            from protrend.pipeline import Pipeline
            Pipeline.register_transformer(cls, **kwargs)

    def __init__(self,
                 transform_stack: Dict[str, str] = None,
                 source: str = None,
                 version: str = None,
                 node: Type[Node] = None,
                 order: int = None):

        """
        The transform object must implement several transformation procedures for a given neomodel Node.
        Thus, the transform must be able to read staging area files and digest/process/transform/integrate/write
        them into fully structured Node instances.

        A TransformerSettings object must be provided to a transformer.
        These settings set up a given data source contained in the staging area for the transformation to take place.

        Pandas is the main engine to read, process, transform and integrate data
        contained in the staging area files into structured nodes.

        :type transform_stack: Dict[str, str]
        :type source: str
        :type version: str
        :type node: Type[Node]
        :type order: int

        :param transform_stack: Dictionary containing the pair name and file name.
        The key should be used to identify the file in the transform stack,
        whereas the value should be the file name in the staging area or data lake
        :param source: The name of the data source in the staging area (e.g. regprecise, collectf, etc)
        :param version: The version of the data source in the staging area (e.g. 0.0.0, 0.0.1, etc)
        :param node: The node type associated with this transformer.
        Note that, it should be created only one transformer for each node
        :param order: The order value is used by the director to rearrange transformers and their execution.
        Transformers having higher order are executed first.
        """

        self.source = source
        self.version = version
        self.node = node
        self.order = order

        if not transform_stack:
            transform_stack = self.default_transform_stack

        self._transform_stack = self.build_stack(transform_stack)

        self._write_stack = []

    def build_stack(self, stack_to_load: Dict[str, str]) -> Dict[str, str]:

        loaded_stack = {}

        for key, file in stack_to_load.items():

            sa_file = os.path.join(Settings.staging_area, self.source, self.version, file)

            if os.path.exists(sa_file):
                loaded_stack[key] = sa_file

            else:
                dl_file = os.path.join(Settings.data_lake, self.source, self.version, file)

                loaded_stack[key] = dl_file

        return loaded_stack

    # --------------------------------------------------------
    # Static properties
    # --------------------------------------------------------
    @property
    def node_factors(self) -> Dict[str, Callable]:
        return self.node.node_factors

    @property
    def node_factors_keys(self) -> List[str]:
        return list(self.node.node_factors.keys())

    @property
    def transform_stack(self) -> Dict[str, str]:
        return self._transform_stack

    @property
    def write_path(self) -> str:
        return os.path.join(Settings.data_lake, self.source, self.version)

    # --------------------------------------------------------
    # Python API
    # --------------------------------------------------------
    def __str__(self):
        return f'{self.__class__.__name__}: {self.source} - {self.version} - ' \
               f'{self.node.node_name()} - {self.node_factors_keys}'

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other: 'Transformer'):

        other_values = {other.__class__.__name__, other.source, other.version, other.node.node_name(),
                        other.node_factors_keys}

        if self.__class__.__name__ not in other_values:
            return False

        if self.source not in other_values:
            return False

        if self.version not in other_values:
            return False

        if self.node.node_name() not in other_values:
            return False

        if self.node_factors not in other_values:
            return False

        return True

    # --------------------------------------------------------
    # Transformer API
    # --------------------------------------------------------
    @abstractmethod
    def transform(self):

        """
        The method responsible for transforming multiple pandas DataFrames into an annotated, cleaned and standardized
        pandas DataFrame ready to be integrated and transformed into structured nodes.

        Interface implementation with unknown signature.
        Concrete implementations are available at the transformer children.

        :return:
        """

        pass

    def integrate(self, df: pd.DataFrame):

        """
        The method responsible for integrating the transformed pandas DataFrame with the current state of the neo4j
        database using neomodel and node class extension

        The node factors are node attributes used for data normalization and integration.
        The node factors are used to retrieve the protrend identifier associated with each row of the database.
        Pandas boolean masks are used to query the database improving search performance.

        The nodes identified by the distinct node factors are marked as update nodes. Thus, these nodes identified by
        their protrend identifier will be updated with new data.

        The nodes missing a hit with the current state of the database will be created according
        to the data contained in the DataFrame.

        :param df: transformed pandas DataFrame
        :return: it creates a new pandas DataFrame of the integrated data
        """
        # take a db snapshot for the current node and ensure uniqueness
        view = self.node_view()
        view = self.standardize_factors(df=view)
        view = self.drop_empty_string(view, *self.node_factors_keys)
        view = self.drop_duplicates(df=view, subset=self.node_factors_keys)

        # ensure uniqueness
        df = self.standardize_factors(df=df)
        df = self.drop_empty_string(df, *self.node_factors_keys)
        df = self.drop_duplicates(df=df, subset=self.node_factors_keys)

        # assign the integration and load columns
        df = df.assign(protrend_id=None, load=None, what='nodes')

        # try to integrate the new nodes by the node factors
        for factor in self.node_factors_keys:
            mapper = view.dropna(subset=[factor])
            mapper = mapper.set_index(factor)
            mapper = mapper['protrend_id']
            to_fill = df[factor].map(mapper)
            filled = df['protrend_id'].fillna(to_fill)
            df = df.assign(protrend_id=filled)

        # those that mismatch the integration by the node factors (that is, the protrend id is missing) will be assigned
        # for creation
        mask = df['protrend_id'].isna()
        df.loc[mask, 'load'] = 'create'
        df.loc[~mask, 'load'] = 'update'

        # inferring the batch size, creating the new ids and assign to the df
        nodes = view['protrend_id'].to_list()
        last_node_idx = self.last_node_index(nodes)
        batch_size = mask.sum()

        batch_ids = self.protrend_identifiers_batch(last_node_idx=last_node_idx, size=batch_size)
        df.loc[mask, 'protrend_id'] = batch_ids

        self.stack_integrated_nodes(df)
        self.stack_nodes(df)

    def write(self):

        if not os.path.exists(self.write_path):
            os.makedirs(self.write_path)

        for json_partial in self._write_stack:
            json_partial()

        self._write_stack = []

    # ----------------------------------------
    # Annotations
    # ----------------------------------------
    @staticmethod
    def annotate_genes(df: pd.DataFrame) -> pd.DataFrame:

        input_values = get_values(df, 'input_values')
        loci = get_values(df, 'locus_tag')
        names = get_values(df, 'name')
        taxa = get_values(df, 'taxonomy')
        uniprot_proteins = get_values(df, 'uniprot_accession')
        ncbi_proteins = get_values(df, 'ncbi_protein')
        ncbi_genbanks = get_values(df, 'genbank_accession')
        ncbi_refseqs = get_values(df, 'refseq_accession')
        ncbi_genes = get_values(df, 'ncbi_gene')

        genes = [GeneDTO(input_value=input_value) for input_value in input_values]

        annotate_genes(dtos=genes,
                       loci=loci,
                       names=names,
                       taxa=taxa,
                       uniprot_proteins=uniprot_proteins,
                       ncbi_proteins=ncbi_proteins,
                       ncbi_genbanks=ncbi_genbanks,
                       ncbi_refseqs=ncbi_refseqs,
                       ncbi_genes=ncbi_genes)

        genes_dict = [dto.to_dict() for dto in genes]

        genes_df = pd.DataFrame(genes_dict)

        strand_mask = (genes_df['strand'] != 'reverse') & (genes_df['strand'] != 'forward')
        genes_df.loc[strand_mask, 'strand'] = None

        return genes_df

    @staticmethod
    def annotate_organisms(df: pd.DataFrame) -> pd.DataFrame:

        input_values = get_values(df, 'input_values')
        identifiers = get_values(df, 'ncbi_taxonomy')
        names = get_values(df, 'name')

        organisms = [OrganismDTO(input_value=input_value) for input_value in input_values]

        annotate_organisms(dtos=organisms, identifiers=identifiers, names=names)

        organisms_dict = [dto.to_dict() for dto in organisms]

        organisms_df = pd.DataFrame(organisms_dict)
        return organisms_df

    # ----------------------------------------
    # Utilities
    # ----------------------------------------
    @classmethod
    def infer_write_stack(cls) -> WriteStack:
        node_name = cls.node.get_default(cls).node_name()
        write_stack = WriteStack(transformed=f'transformed_{node_name}',
                                 integrated=f'integrated_{node_name}',
                                 nodes=f'nodes_{node_name}',
                                 connected=None)
        return write_stack

    def empty_frame(self) -> pd.DataFrame:
        cols = [col for col in self.columns if col != 'protrend_id']
        return pd.DataFrame(columns=cols)

    def integrated_empty_frame(self) -> pd.DataFrame:
        cols = [col for col in self.columns]
        cols += ['load', 'what', 'node']
        return pd.DataFrame(columns=cols)

    def stack_json(self, name: str, df: pd.DataFrame):
        df = df.copy()
        df = df.reset_index(drop=True)
        fp = os.path.join(self.write_path, f'{name}.json')
        json_partial = partial(write_json_frame, file_path=fp, df=df)
        self._write_stack.append(json_partial)

    def stack_nodes(self, df: pd.DataFrame):
        if df.empty:
            df = self.integrated_empty_frame()

        else:
            df.loc[:, 'node'] = self.node.node_name()

        node_cols = self.node_factors_keys + ['load', 'what', 'node']
        cols_to_drop = [col for col in df.columns if col not in node_cols]
        df = df.drop(columns=cols_to_drop)

        df_name = f'nodes_{self.node.node_name()}'
        self.stack_json(df_name, df)

    def stack_integrated_nodes(self, df: pd.DataFrame):
        if df.empty:
            df = self.integrated_empty_frame()

        else:
            df.loc[:, 'node'] = self.node.node_name()

        df_name = f'integrated_{self.node.node_name()}'
        self.stack_json(df_name, df)

    def stack_transformed_nodes(self, df: pd.DataFrame):
        if df.empty:
            df = self.empty_frame()
        df_name = f'transformed_{self.node.node_name()}'
        self.stack_json(df_name, df)

    @staticmethod
    def contact_stacks(stack: Dict[str, str], taxa: Dict[str, str], **kwargs) -> pd.DataFrame:
        dfs = []
        for file in stack:
            df = read_from_stack(stack=stack,
                                 file=file,
                                 **kwargs)
            df['taxonomy'] = taxa[file]
            dfs.append(df)

        return pd.concat(dfs)

    def standardize_factors(self, df: pd.DataFrame) -> pd.DataFrame:
        df = apply_processors(df, **self.node_factors)
        return df

    @staticmethod
    def drop_empty_string(df: pd.DataFrame, *cols: str) -> pd.DataFrame:
        n_rows, _ = df.shape

        mask = pd.Series([False] * n_rows)

        cols_mask = [df[col] != '' for col in cols]

        if cols_mask:
            mask = pd.concat(cols_mask, axis=1)
            mask = mask.all(axis=1)

        new_df = df[mask].copy()
        return new_df

    @staticmethod
    def create_input_value(df: pd.DataFrame, col: str) -> pd.DataFrame:
        df['input_value'] = df[col]
        return df

    @staticmethod
    def select_columns(df: pd.DataFrame, *columns: str) -> pd.DataFrame:
        df = df.copy()
        df = df.loc[:, list(columns)]
        return df

    @staticmethod
    def merge_columns(df: pd.DataFrame, column: str, left: str, right: str) -> pd.DataFrame:
        df[column] = df[left].fillna(df[right])
        df = df.drop(columns=[left, right])
        return df

    @staticmethod
    def concat_columns(df: pd.DataFrame, column: str, left: str, right: str) -> pd.DataFrame:
        df[column] = df[left] + df[right]
        df = df.drop(columns=[left, right])
        return df

    @staticmethod
    def group_by(df: pd.DataFrame,
                 column: str,
                 aggregation: Dict[str, Callable],
                 default: Callable = take_last) -> pd.DataFrame:

        agg = {col: default for col in df.columns if col != column}
        agg.update(aggregation)

        df = df.groupby(df[column]).aggregate(agg)
        df = df.reset_index()
        return df

    @staticmethod
    def drop_duplicates(df: pd.DataFrame,
                        subset: List[str],
                        perfect_match: bool = False,
                        preserve_nan: bool = True) -> pd.DataFrame:

        if perfect_match and preserve_nan:

            df = df[(~df.duplicated(subset=subset)) | df[subset].isnull().any(axis=1)]

        elif perfect_match and not preserve_nan:

            df = df.drop_duplicates(subset=subset)

        elif not perfect_match and preserve_nan:

            for col in subset:
                df = df[(~df.duplicated(subset=[col])) | df[col].isnull()]

        elif not perfect_match and not preserve_nan:

            for col in subset:
                df = df.drop_duplicates(subset=[col])

        df = df.reset_index(drop=True)

        return df

    def interaction_hash(self, df: pd.DataFrame) -> pd.DataFrame:
        # filter by organism + regulator + gene + tfbs + effector + regulatory effect
        df2 = apply_processors(df,
                               organism=to_list_nan,
                               regulator=to_list_nan,
                               gene=to_list_nan,
                               tfbs=to_list_nan,
                               effector=to_list_nan,
                               regulatory_effect=to_list_nan)

        ri_series_hash = df2['organism'].copy()
        ri_series_hash += df2['regulator'].copy()
        ri_series_hash += df2['gene'].copy()
        ri_series_hash += df2['tfbs'].copy()
        ri_series_hash += df2['effector'].copy()
        ri_series_hash += df2['regulatory_effect'].copy()

        df = df.assign(regulatory_interaction_hash=ri_series_hash)
        df = apply_processors(df, regulatory_interaction_hash=regulatory_interaction_hash)
        df = self.drop_duplicates(df=df, subset=['regulatory_interaction_hash'], perfect_match=True)
        df = df.dropna(subset=['regulatory_interaction_hash'])

        return df

    def protrend_identifiers_batch(self, last_node_idx: int, size: int) -> List[str]:

        return [protrend_id_encoder(self.node.header, self.node.entity, i)
                for i in range(last_node_idx + 1, last_node_idx + size + 1)]

    @staticmethod
    def last_node_index(nodes: List[str]) -> int:
        if not nodes:
            return 0

        sorted_nodes = sorted(nodes, key=protrend_id_decoder, reverse=True)
        return protrend_id_decoder(sorted_nodes[0])

    def node_view(self) -> pd.DataFrame:
        df = self.node.node_to_df()
        return df
