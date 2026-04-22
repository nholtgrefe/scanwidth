"""Shared utility helpers for scanwidth algorithms."""

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
    """Return a sentinel value larger than any possible edge- or node-scanwidth.

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


def sblocks(graph: nx.DiGraph) -> List[Set]:
    """Return node sets for s-blocks in merge order.

    Parameters
    ----------
    graph : nx.DiGraph
        Input graph.

    Returns
    -------
    List[Set]
        Ordered list of s-block node sets.
    """
    roots = {v for v in graph.nodes() if graph.in_degree(v) == 0}
    aux = graph.to_undirected()

    for root1 in roots:
        for root2 in roots:
            if (
                root1 != root2
                and (root1, root2) not in aux.edges()
                and (root2, root1) not in aux.edges()
            ):
                aux.add_edge(root1, root2)

    sblock_sets = list(nx.biconnected_components(aux))
    dcut_vertices = list(nx.articulation_points(aux))

    sblock_cut_tree = nx.Graph()
    for v in dcut_vertices:
        sblock_cut_tree.add_node(v)

    rootblock_index = None
    for i, block in enumerate(sblock_sets):
        node_name = f"block_{i}"
        if roots.issubset(block):
            rootblock_index = i
        sblock_cut_tree.add_node(node_name)
        for v in dcut_vertices:
            if v in block:
                sblock_cut_tree.add_edge(v, node_name)

    if rootblock_index is None:
        return sblock_sets

    sblock_order = list(
        nx.dfs_preorder_nodes(
            sblock_cut_tree, source=f"block_{rootblock_index}",
        )
    )
    sblock_order = [
        int(name[6:])
        for name in sblock_order
        if str(name).startswith("block_")
    ]
    return [sblock_sets[i] for i in sblock_order]


def chain_vertices(graph: nx.DiGraph) -> List:
    """Return degree-(1,1) chain vertices in ``graph``.

    Notes
    -----
    Chain vertices are also called flow vertices.

    Parameters
    ----------
    graph : nx.DiGraph
        Input graph.

    Returns
    -------
    List
        Vertices with indegree 1 and outdegree 1.
    """
    return [
        v for v in graph.nodes()
        if graph.out_degree(v) == 1 and graph.in_degree(v) == 1
    ]


def is_directed_tree(graph: nx.DiGraph) -> bool:
    """Return whether ``graph`` is a directed rooted tree.

    Parameters
    ----------
    graph : nx.DiGraph
        Input graph.

    Returns
    -------
    bool
        True if the graph is a directed rooted tree, else False.
    """
    n = graph.number_of_nodes()
    if n == 0:
        return True
    if graph.number_of_edges() != n - 1:
        return False
    if not nx.is_weakly_connected(graph):
        return False
    if not nx.is_directed_acyclic_graph(graph):
        return False
    roots = [v for v in graph.nodes() if graph.in_degree(v) == 0]
    return len(roots) == 1

