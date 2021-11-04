from typing import List, Any, Set, Generator, Tuple

from networkx import Graph, connected_components


def build_graph(nodes: List[List[Any]]) -> Graph:
    graph_net = Graph()

    for nodes_list in nodes:
        # each nodes_list is a bunch of nodes
        graph_net.add_nodes_from(nodes_list)

        # associated with a set of edges
        edges = build_edges(nodes_list)
        graph_net.add_edges_from(edges)

    return graph_net


def build_edges(nodes: List[Any]) -> Generator[Tuple[Any, Any], Any, Any]:
    forward = [None] + nodes
    backwards = nodes + [None]

    for current, last in zip(forward, backwards):

        if current is None or last is None:
            continue

        yield current, last


def find_connected_nodes(graph_net: Graph) -> List[Set[Any]]:

    return list(connected_components(graph_net))
