"""Shared internal helpers for edge-scanwidth solvers."""

from __future__ import annotations

from typing import List, Set, Union

import networkx as nx


def delta_in(
    graph: nx.DiGraph,
    vertex_set: Union[Set, List],
    sink: bool = True,
) -> int:
    """Return the indegree of ``vertex_set`` in ``graph``.

    Parameters
    ----------
    graph : nx.DiGraph
        Directed graph containing ``vertex_set``.
    vertex_set : Union[Set, List]
        Vertices for which to compute the indegree.
    sink : bool, optional
        If True, treat ``vertex_set`` as a sinkset to speed up computation.

    Returns
    -------
    int
        The indegree of ``vertex_set``.
    """
    res = 0
    if sink:
        if len(vertex_set) < len(graph.nodes()) / 2:
            for v in vertex_set:
                res = res + graph.in_degree(v) - graph.out_degree(v)
        else:
            for v in graph.nodes:
                if v not in vertex_set:
                    res = res - graph.in_degree(v) + graph.out_degree(v)
        return res

    for (u, v) in graph.edges():
        if u not in vertex_set and v in vertex_set:
            res = res + 1
    return res


def find_component(components: List[Set], v: object) -> Set:
    """Return the component in ``components`` that contains ``v``.

    Parameters
    ----------
    components : List[Set]
        List of disjoint vertex sets.
    v : object
        Vertex to look up.

    Returns
    -------
    Set
        The component containing ``v``, or an empty set if none found.
    """
    for comp in components:
        if v in comp:
            return comp
    return set()


def infinity_for(graph: nx.DiGraph) -> int:
    """Return a sentinel value larger than any possible edge scanwidth.

    Parameters
    ----------
    graph : nx.DiGraph
        Input graph.

    Returns
    -------
    int
        Sentinel infinity for the graph.
    """
    return graph.number_of_edges() + 1
