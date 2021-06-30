from typing import List, Type

import pandas as pd

from protrend.models.node import Node


class NodeImporter:

    def __init__(self, node_cls: Type[Node], path: str):
        self.node_cls: Type[Node] = node_cls
        self.path: str = path
        self._nodes: list = []

    @property
    def nodes(self) -> List[Node]:
        return self._nodes

    def append(self, node: Node):
        self._nodes.append(node)

    @property
    def id(self) -> str:
        return self.node_cls.property_as_id

    @property
    def name(self) -> str:
        return self.node_cls.cls_name()

    @property
    def id_column(self) -> str:
        return f'{self.id}:ID({self.name})'

    @property
    def node_columns(self) -> List[str]:
        return self.node_cls.cls_keys()

    @property
    def node_csv(self) -> str:
        return fr'{self.path}\{self.name}.csv'

    def relationship_csv(self, relationship) -> str:
        rel_type = relationship.definition['relation_type']
        out_name = relationship.definition['node_class'].cls_name()
        return fr'{self.path}\{self.name}_{rel_type}_{out_name}.csv'

    def _build_node_import(self) -> pd.DataFrame:
        series = [node.to_series() for node in self._nodes]
        df = pd.DataFrame(series, columns=self.node_columns)
        df.rename(columns={self.id: self.id_column}, inplace=True)
        return df

    def _build_relationship_import(self, name, relationship) -> pd.DataFrame:
        dfs = [node.relationship_to_df(name) for node in self.nodes]

        if dfs:
            df = pd.concat(dfs)

        else:
            df = pd.DataFrame(columns=['start', 'end'])

        out_name = relationship.definition['node_class'].cls_name()
        df.rename(columns={'start': f':START_ID({self.name})', 'end': f':END_ID({out_name})'},
                  inplace=True)
        return df

    def build_imports(self):

        if not self.nodes:
            return

        df = self._build_node_import()
        df.to_csv(self.node_csv, index=False)

        for key, val in self.node_cls.cls_relationships().items():
            df = self._build_relationship_import(key, val)
            csv_file = self.relationship_csv(val)
            df.to_csv(csv_file, index=False)

    @property
    def args(self) -> List[str]:

        if not self.nodes:
            return []

        args = [f'--nodes={self.name}={self.node_csv}']

        for key, val in self.node_cls.cls_relationships().items():
            csv_file = self.relationship_csv(val)
            rel_type = val.definition['relation_type']
            arg = f'--relationships={rel_type}={csv_file}'
            args.append(arg)

        return args
