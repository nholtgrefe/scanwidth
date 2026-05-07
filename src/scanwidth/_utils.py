"""Shared utility helpers for scanwidth algorithms."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Set, Tuple, Union
from weakref import WeakKeyDictionary

import networkx as nx


@dataclass
class _ParentBoundaryCache:
    """Per-graph cache for fast distinct-parent boundary counting."""

    node_to_idx: dict
    preds_idx: List[List[int]]
    succ_idx: List[List[int]]
    in_marks: List[int]
    parent_marks: List[int]
    in_epoch: int = 1
    parent_epoch: int = 1


_PARENT_BOUNDARY_CACHE: "WeakKeyDictionary[nx.DiGraph, _ParentBoundaryCache]" = (
    WeakKeyDictionary()
)

_WCC_PARTITION_CACHE: "WeakKeyDictionary[nx.DiGraph, dict[frozenset, Tuple[frozenset, ...]]]" = (
    WeakKeyDictionary()
)


def _parent_boundary_cache(graph: nx.DiGraph) -> _ParentBoundaryCache:
    """Return cached integer adjacency structures for ``graph``."""
    cache = _PARENT_BOUNDARY_CACHE.get(graph)
    if cache is not None:
        return cache

    nodes = list(graph.nodes())
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    n = len(nodes)

    preds_idx: List[List[int]] = [[] for _ in range(n)]
    succ_idx: List[List[int]] = [[] for _ in range(n)]
    for (u, v) in graph.edges():
        ui = node_to_idx[u]
        vi = node_to_idx[v]
        succ_idx[ui].append(vi)
        preds_idx[vi].append(ui)

    cache = _ParentBoundaryCache(
        node_to_idx=node_to_idx,
        preds_idx=preds_idx,
        succ_idx=succ_idx,
        in_marks=[0] * n,
        parent_marks=[0] * n,
    )
    _PARENT_BOUNDARY_CACHE[graph] = cache
    return cache


def _next_epoch(current: int) -> int:
    """Return next positive epoch value."""
    next_value = current + 1
    if next_value <= 0:
        return 1
    return next_value


def delta_in(
    graph: nx.DiGraph,
    vertex_set: Union[Set, List],
    sink: bool = False,
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


def induced_subgraph_roots(
    graph: nx.DiGraph,
    vertices: Union[Set, List],
) -> List:
    """Return roots of the induced subgraph ``G[W]`` for ``W = vertices``.

    A vertex ``v \\in W`` is a root of ``G[W]`` iff ``graph`` has no arc
    ``(u, v)`` with ``u \\in W``.

    Parameters
    ----------
    graph : nx.DiGraph
        Host graph.
    vertices : Union[Set, List]
        Vertex set ``W``. If a list, roots are returned in that list order;
        if a set, order follows ``list(vertices)``.

    Returns
    -------
    List
        Roots of ``G[W]``, in the order induced by ``list(vertices)``.
    """
    vertex_list = list(vertices)
    W = set(vertex_list)
    return [
        v for v in vertex_list
        if not any(pred in W for pred in graph.predecessors(v))
    ]


def induced_weakly_connected_components(
    graph: nx.DiGraph,
    vertices: Union[Set, List],
) -> List[Set]:
    """Return weakly connected components of the induced subgraph ``G[W]``.

    Components are returned as mutable ``set`` objects (like NetworkX).
    Results are cached per ``(graph, frozenset(W))`` so repeated queries
    for the same vertex set reuse one NetworkX decomposition.

    Component lists are sorted deterministically (lexicographic order of
    sorted string vertex labels within each component) so cache hits and
    misses agree on iteration order.

    Parameters
    ----------
    graph : nx.DiGraph
        Host graph.
    vertices : Union[Set, List]
        Vertex set ``W`` inducing ``G[W]``.

    Returns
    -------
    List[Set]
        Weakly connected components of ``G[W]``.
    """
    W = frozenset(vertices)
    if not W:
        return []

    inner = _WCC_PARTITION_CACHE.get(graph)
    if inner is None:
        inner = {}
        _WCC_PARTITION_CACHE[graph] = inner

    cached = inner.get(W)
    if cached is not None:
        return [set(c) for c in cached]

    subgraph = graph.subgraph(W)
    raw = list(nx.weakly_connected_components(subgraph))
    parts = tuple(
        frozenset(c)
        for c in sorted(
            raw,
            key=lambda s: tuple(sorted(str(v) for v in s)),
        )
    )
    inner[W] = parts
    return [set(c) for c in parts]


def delta_in_parents(
    graph: nx.DiGraph,
    vertex_set: Union[Set, List],
    sink: bool = False,
) -> int:
    """Return number of distinct parents outside ``vertex_set`` with an edge into it.

    Parameters
    ----------
    graph : nx.DiGraph
        Directed graph containing ``vertex_set``.
    vertex_set : Union[Set, List]
        Vertices that define the sink-side set.
    sink : bool, optional
        If True, treat ``vertex_set`` as a sink-set and use an optimized
        strategy that inspects only one side of the boundary. If False,
        scan all edges directly.

    Returns
    -------
    int
        Number of distinct parent vertices outside ``vertex_set`` that
        have at least one edge into ``vertex_set``.
    """
    cache = _parent_boundary_cache(graph)
    node_to_idx = cache.node_to_idx
    n = len(cache.in_marks)

    in_epoch = _next_epoch(cache.in_epoch)
    cache.in_epoch = in_epoch
    in_marks = cache.in_marks
    vertices_idx: List[int] = []
    for vertex in vertex_set:
        idx = node_to_idx[vertex]
        vertices_idx.append(idx)
        in_marks[idx] = in_epoch

    parent_epoch = _next_epoch(cache.parent_epoch)
    cache.parent_epoch = parent_epoch
    parent_marks = cache.parent_marks

    boundary_count = 0
    if sink:
        # Pick the smaller side of the cut for boundary detection.
        if len(vertices_idx) * 2 < n:
            # ``W`` smaller: aggregate outside parents from incoming arcs of ``W``.
            for vertex_idx in vertices_idx:
                for parent_idx in cache.preds_idx[vertex_idx]:
                    if (
                        in_marks[parent_idx] != in_epoch
                        and parent_marks[parent_idx] != parent_epoch
                    ):
                        parent_marks[parent_idx] = parent_epoch
                        boundary_count += 1
            return boundary_count

        # ``V \ W`` smaller: count outside vertices with at least one child in ``W``.
        for parent_idx in range(n):
            if in_marks[parent_idx] == in_epoch:
                continue
            for child_idx in cache.succ_idx[parent_idx]:
                if in_marks[child_idx] == in_epoch:
                    boundary_count += 1
                    break
        return boundary_count

    for parent_idx, child_idx in (
        (node_to_idx[u], node_to_idx[v]) for (u, v) in graph.edges()
    ):
        if (
            in_marks[parent_idx] != in_epoch
            and in_marks[child_idx] == in_epoch
            and parent_marks[parent_idx] != parent_epoch
        ):
            parent_marks[parent_idx] = parent_epoch
            boundary_count += 1
    return boundary_count


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


def sblocks(graph: nx.DiGraph) -> List[Tuple[Set, bool]]:
    """Return s-blocks as ``(block, is_root_block)`` in merge order.

    Parameters
    ----------
    graph : nx.DiGraph
        Input graph.

    Notes
    -----
    If the graph has multiple roots, the decomposition contains a special
    root block. This block is not biconnected in the original graph and
    can be identified by its corresponding ``True`` flag.

    Returns
    -------
    List[Tuple[Set, bool]]
        Ordered list of tuples ``(block_nodes, is_root_block)``.
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
        block_node = ("block", i)
        if roots.issubset(block):
            rootblock_index = i
        sblock_cut_tree.add_node(block_node)
        for v in dcut_vertices:
            if v in block:
                sblock_cut_tree.add_edge(v, block_node)

    if rootblock_index is None:
        return [(block, False) for block in sblock_sets]

    sblock_order = list(
        nx.dfs_preorder_nodes(
            sblock_cut_tree, source=("block", rootblock_index),
        )
    )
    sblock_order = [
        node[1]
        for node in sblock_order
        if isinstance(node, tuple) and len(node) == 2 and node[0] == "block"
    ]
    has_multiple_roots = len(roots) > 1
    return [
        (sblock_sets[i], has_multiple_roots and i == rootblock_index)
        for i in sblock_order
    ]


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

