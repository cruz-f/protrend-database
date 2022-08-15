from typing import Dict, Type

import pandas as pd

from protrend.model import *


entities = {
    'effector': Effector,
    'evidence': Evidence,
    'gene': Gene,
    'motif': Motif,
    'operon': Operon,
    'organism': Organism,
    'pathway': Pathway,
    'publication': Publication,
    'regulator': Regulator,
    'regulatory family': RegulatoryFamily,
    'regulatory interaction': RegulatoryInteraction,
    'source': Source,
    'tfbs': TFBS,
}


def database_short_report(nodes: Dict[str, Type[BaseNode]] = None) -> pd.DataFrame:
    """
    Generate a report of the database.
    """
    if not nodes:
        nodes = entities.copy()

    results = {}
    for name, node in nodes.items():
        results[name] = [len(node.nodes), len(node.node_properties()), len(node.node_relationships())]

    return pd.DataFrame(results)


def database_long_report(nodes: Dict[str, Type[BaseNode]] = None) -> pd.DataFrame:
    """
    Generate a long report of the database.
    """
    if not nodes:
        nodes = entities.copy()

    results = []
    keys = []
    for name, node in nodes.items():
        df = node.node_to_df()
        results.append(df)
        keys.append(name)

    return pd.concat(results, keys=keys, axis=0)


if __name__ == '__main__':
    short_report = database_short_report()
    long_report = database_long_report()
