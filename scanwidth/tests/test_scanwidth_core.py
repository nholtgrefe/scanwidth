"""Core correctness tests for scanwidth computations.

The tests in this module focus on small hand-crafted DAGs where expected
scanwidth values are easy to verify analytically.
"""

from __future__ import annotations

import networkx as nx

from scanwidth import DAG, Extension, edge_scanwidth, TreeExtension


def _build_chain() -> nx.DiGraph:
    """Build a simple 3-node chain DAG.

    Returns
    -------
    nx.DiGraph
        Directed graph with edges (1, 2), (2, 3).
    """
    return nx.DiGraph([(1, 2), (2, 3)])


def _build_two_to_one() -> nx.DiGraph:
    """Build a two-sources-to-one-sink DAG.

    Returns
    -------
    nx.DiGraph
        Directed graph with edges (1, 3), (2, 3).
    """
    return nx.DiGraph([(1, 3), (2, 3)])


def _build_star_to_sink(k: int) -> nx.DiGraph:
    """Build a DAG with ``k`` source vertices pointing to one sink.

    Parameters
    ----------
    k : int
        Number of source vertices.

    Returns
    -------
    nx.DiGraph
        Directed graph with edges ``(1, sink), ..., (k, sink)`` where
        ``sink = k + 1``.
    """
    sink = k + 1
    edges = [(i, sink) for i in range(1, k + 1)]
    return nx.DiGraph(edges)


def test_extension_scanwidth_known_value() -> None:
    """Extension scanwidth matches known value on a small DAG."""
    graph = _build_two_to_one()
    extension = Extension(dag=DAG(graph), ordering=[3, 1, 2])

    assert len(extension.edge_scanwidth_bag(3)) == 2
    assert extension.node_scanwidth_bag(3) == {1, 2}
    assert extension.edge_scanwidth() == 2
    assert extension.node_scanwidth() == 2


def test_tree_extension_roundtrip_preserves_scanwidth() -> None:
    """Canonical tree extension preserves scanwidth and validity."""
    graph = _build_two_to_one()
    extension = Extension(dag=DAG(graph), ordering=[3, 1, 2])

    tree_extension = extension.canonical_tree_extension()
    roundtrip_extension = tree_extension.to_extension()

    assert tree_extension.is_canonical()
    assert tree_extension.edge_scanwidth() == extension.edge_scanwidth()
    assert tree_extension.node_scanwidth() == extension.node_scanwidth()
    assert roundtrip_extension.edge_scanwidth() == extension.edge_scanwidth()
    assert roundtrip_extension.node_scanwidth() == extension.node_scanwidth()


def test_extension_properties_and_canonical_alias() -> None:
    """Extension properties and canonical alias behave as expected."""
    dag = DAG(_build_two_to_one())
    extension = Extension(dag=dag, ordering=[3, 1, 2])

    # dag property exposes the wrapped DAG object
    assert extension.dag is dag

    # ordering property returns a copy
    ordering_view = extension.ordering
    ordering_view.append("x")
    assert extension.ordering == [3, 1, 2]

    # both canonical constructors are equivalent interfaces
    tree_via_to = extension.to_canonical_tree_extension()
    tree_via_alias = extension.canonical_tree_extension()
    assert tree_via_to.is_canonical()
    assert tree_via_alias.is_canonical()
    assert tree_via_to.edge_scanwidth() == tree_via_alias.edge_scanwidth()


def test_canonical_tree_extension_preserves_vertex_bags() -> None:
    """Canonical tree extension matches extension bags per vertex."""
    extension = Extension(dag=DAG(_build_two_to_one()), ordering=[3, 1, 2])
    tree_extension = extension.canonical_tree_extension()

    for vertex in extension.ordering:
        assert extension.edge_scanwidth_bag(vertex) == tree_extension.edge_scanwidth_bag(vertex)
        assert extension.node_scanwidth_bag(vertex) == tree_extension.node_scanwidth_bag(vertex)


def test_tree_extension_edge_scanwidth_bag_properties() -> None:
    """Each edge bag contains exactly incoming boundary edges."""
    graph = _build_two_to_one()
    tree_extension = Extension(
        dag=DAG(graph), ordering=[3, 1, 2],
    ).canonical_tree_extension()

    for vertex in tree_extension.tree.nodes():
        bag = tree_extension.edge_scanwidth_bag(vertex)
        left = nx.descendants(tree_extension.tree, vertex) | {vertex}
        expected = {
            (u, w) for (u, w) in tree_extension.dag.graph.edges()
            if u not in left and w in left
        }
        assert bag == expected
        assert len(bag) <= tree_extension.edge_scanwidth()
        assert tree_extension.node_scanwidth_bag(vertex) == {u for (u, _) in bag}


def test_tree_extension_properties_and_roundtrip_ordering() -> None:
    """TreeExtension properties and roundtrip ordering are consistent."""
    extension = Extension(dag=DAG(_build_two_to_one()), ordering=[3, 1, 2])
    tree_extension = extension.canonical_tree_extension()
    roundtrip = tree_extension.to_extension()

    assert tree_extension.dag is extension.dag
    assert set(tree_extension.tree.nodes()) == set(extension.ordering)
    assert roundtrip.ordering == extension.ordering


def test_tree_extension_edge_scanwidth_bag_rejects_unknown_vertex() -> None:
    """Bag query fails for vertices outside the tree extension."""
    graph = _build_two_to_one()
    tree_extension = Extension(
        dag=DAG(graph), ordering=[3, 1, 2],
    ).canonical_tree_extension()

    try:
        _ = tree_extension.edge_scanwidth_bag("missing")
        assert False, "Expected ValueError for unknown tree vertex."
    except ValueError:
        pass


def test_extension_edge_scanwidth_bag_rejects_unknown_vertex() -> None:
    """Bag query fails for vertices outside the extension order."""
    extension = Extension(dag=DAG(_build_two_to_one()), ordering=[3, 1, 2])

    try:
        _ = extension.edge_scanwidth_bag("missing")
        assert False, "Expected ValueError for unknown extension vertex."
    except ValueError:
        pass


def test_extension_node_scanwidth_bag_matches_edge_parents() -> None:
    """Node bag equals parent set of the edge bag per extension vertex."""
    extension = Extension(dag=DAG(_build_two_to_one()), ordering=[3, 1, 2])
    for vertex in extension.ordering:
        edge_bag = extension.edge_scanwidth_bag(vertex)
        assert extension.node_scanwidth_bag(vertex) == {u for (u, _) in edge_bag}
        assert len(extension.node_scanwidth_bag(vertex)) <= extension.node_scanwidth()


def test_tree_extension_node_scanwidth_bag_matches_edge_parents() -> None:
    """Node bag equals parent set of edge bag per tree extension vertex."""
    tree_extension = Extension(
        dag=DAG(_build_two_to_one()), ordering=[3, 1, 2],
    ).canonical_tree_extension()
    for vertex in tree_extension.tree.nodes():
        edge_bag = tree_extension.edge_scanwidth_bag(vertex)
        assert tree_extension.node_scanwidth_bag(vertex) == {u for (u, _) in edge_bag}
        assert len(tree_extension.node_scanwidth_bag(vertex)) <= tree_extension.node_scanwidth()


def test_exact_functions_match_known_chain_scanwidth() -> None:
    """All exact solvers compute scanwidth 1 on a directed chain."""
    dag = DAG(_build_chain())

    sw_exhaustive, ext_exhaustive = edge_scanwidth(
        dag, algorithm="exhaustive", reduce=False,
    )
    sw_two_partition, ext_two_partition = edge_scanwidth(
        dag, algorithm="two_partition", reduce=False,
    )
    sw_partition, ext_partition = edge_scanwidth(
        dag, algorithm="three_partition", reduce=False,
    )
    sw_xp, ext_xp = edge_scanwidth(dag, algorithm="xp", reduce=False)

    assert sw_exhaustive == sw_two_partition == sw_partition == sw_xp == 1
    assert ext_exhaustive.edge_scanwidth() == 1
    assert ext_two_partition.edge_scanwidth() == 1
    assert ext_partition.edge_scanwidth() == 1
    assert ext_xp.edge_scanwidth() == 1


def test_heuristics_return_valid_extensions_and_consistent_scanwidth() -> None:
    """Heuristics return consistent scanwidth values and valid extensions."""
    dag = DAG(_build_chain())

    sw_greedy, ext_greedy = edge_scanwidth(dag, algorithm="greedy", reduce=False)
    sw_cut, ext_cut = edge_scanwidth(dag, algorithm="cut_splitting", reduce=False)
    sw_anneal, ext_anneal = edge_scanwidth(
        dag,
        algorithm="simulated_annealing",
        max_iter=5,
        verbose=False,
        reduce=False,
        seed=7,
    )

    assert sw_greedy == ext_greedy.edge_scanwidth()
    assert sw_cut == ext_cut.edge_scanwidth()
    assert sw_anneal == ext_anneal.edge_scanwidth()


def _assert_star_to_sink_scanwidth(expected_scanwidth: int) -> None:
    """Exact algorithm returns known scanwidths for simple star DAGs.

    Parameters
    ----------
    expected_scanwidth : int
        Expected optimal scanwidth of the generated DAG.
    """
    dag = DAG(_build_star_to_sink(expected_scanwidth))
    sw, ext = edge_scanwidth(dag, algorithm="exhaustive", reduce=False)

    assert ext is not None
    assert sw == expected_scanwidth
    assert ext.edge_scanwidth() == expected_scanwidth


def test_star_to_sink_scanwidth_3() -> None:
    """Exact algorithm returns scanwidth 3 for the 3-star DAG."""
    _assert_star_to_sink_scanwidth(3)


def test_star_to_sink_scanwidth_4() -> None:
    """Exact algorithm returns scanwidth 4 for the 4-star DAG."""
    _assert_star_to_sink_scanwidth(4)


def test_edge_scanwidth_entrypoint_xp() -> None:
    """Main entry point computes expected value for XP on a chain."""
    dag = DAG(_build_chain())
    sw, ext = edge_scanwidth(dag, algorithm="xp", reduce=False)

    assert sw == 1
    assert ext.edge_scanwidth() == 1


def test_edge_scanwidth_entrypoint_xp_rejects_fixed_k() -> None:
    """Public XP entrypoint rejects fixed-k mode."""
    dag = DAG(_build_chain())
    try:
        _ = edge_scanwidth(dag, algorithm="xp", k=1, reduce=False)
        assert False, "Expected ValueError when passing k to public xp API."
    except ValueError:
        pass


def test_edge_scanwidth_reduce_matches_no_reduce() -> None:
    """Reduced and unreduced XP produce identical scanwidth on a star DAG."""
    dag = DAG(_build_star_to_sink(4))

    sw_reduced, ext_reduced = edge_scanwidth(
        dag, algorithm="xp", reduce=True,
    )
    sw_direct, ext_direct = edge_scanwidth(
        dag, algorithm="xp", reduce=False,
    )

    assert sw_reduced == sw_direct == 4


def test_edge_scanwidth_reduce_greedy_on_chain() -> None:
    """Greedy with s-block reduction returns a valid extension."""
    dag = DAG(_build_chain())
    sw, ext = edge_scanwidth(dag, algorithm="greedy", reduce=True)

    assert sw == ext.edge_scanwidth()


def test_dag_init_rejects_cyclic_graph() -> None:
    """DAG initialization rejects directed cyclic graphs."""
    cyclic = nx.DiGraph([(1, 2), (2, 1)])

    try:
        _ = DAG(cyclic)
        assert False, "Expected ValueError for cyclic directed graph."
    except ValueError:
        pass


def test_extension_init_rejects_invalid_extension() -> None:
    """Extension initialization rejects invalid node ordering."""
    graph = _build_chain()

    try:
        _ = Extension(dag=DAG(graph), ordering=[1, 2, 3])
        assert False, "Expected ValueError for invalid extension ordering."
    except ValueError:
        pass


def test_extension_init_requires_dag_wrapper() -> None:
    """Extension initialization requires a DAG wrapper as graph input."""
    graph = _build_chain()
    try:
        _ = Extension(dag=graph, ordering=[3, 2, 1])  # type: ignore[arg-type]
        assert False, "Expected TypeError when graph is not a DAG instance."
    except TypeError:
        pass


def test_extension_init_requires_ordering_list() -> None:
    """Extension initialization requires ordering to be a list."""
    graph = _build_chain()
    try:
        _ = Extension(dag=DAG(graph), ordering="3 2 1")  # type: ignore[arg-type]
        assert False, "Expected TypeError when ordering is not a list."
    except TypeError:
        pass


def test_tree_extension_init_rejects_invalid_tree() -> None:
    """TreeExtension initialization rejects invalid tree extensions."""
    graph = _build_chain()
    invalid_tree = nx.DiGraph([(1, 2), (1, 3)])

    try:
        _ = TreeExtension(dag=DAG(graph), tree=invalid_tree)
        assert False, "Expected ValueError for invalid tree extension."
    except ValueError:
        pass


def test_tree_extension_init_requires_dag_wrapper() -> None:
    """TreeExtension initialization requires a DAG wrapper as graph input."""
    graph = _build_chain()
    tree = nx.DiGraph([(1, 2), (2, 3)])
    try:
        _ = TreeExtension(dag=graph, tree=tree)  # type: ignore[arg-type]
        assert False, "Expected TypeError when graph is not a DAG instance."
    except TypeError:
        pass
