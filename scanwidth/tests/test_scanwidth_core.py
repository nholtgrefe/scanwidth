"""Core correctness tests for scanwidth computations.

The tests in this module focus on small hand-crafted DAGs where expected
scanwidth values are easy to verify analytically.
"""

from __future__ import annotations

import networkx as nx

from scanwidth import DAG, Extension
from scanwidth.edge_scanwidth import edge_scanwidth
from scanwidth.tree_extension import TreeExtension


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
    extension = Extension(graph=graph, sigma=[3, 1, 2])

    assert extension.edge_scanwidth_at_vertex_i(0) == 2
    assert extension.edge_scanwidth() == 2


def test_tree_extension_roundtrip_preserves_scanwidth() -> None:
    """Canonical tree extension preserves scanwidth and validity."""
    graph = _build_two_to_one()
    extension = Extension(graph=graph, sigma=[3, 1, 2])

    tree_extension = extension.canonical_tree_extension()
    roundtrip_extension = tree_extension.to_extension()

    assert tree_extension.is_canonical()
    assert tree_extension.edge_scanwidth() == extension.edge_scanwidth()
    assert roundtrip_extension.edge_scanwidth() == extension.edge_scanwidth()


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
        _ = Extension(graph=graph, sigma=[1, 2, 3])
        assert False, "Expected ValueError for invalid extension ordering."
    except ValueError:
        pass


def test_tree_extension_init_rejects_invalid_tree() -> None:
    """TreeExtension initialization rejects invalid tree extensions."""
    graph = _build_chain()
    invalid_tree = nx.DiGraph([(1, 2), (1, 3)])

    try:
        _ = TreeExtension(graph, invalid_tree)
        assert False, "Expected ValueError for invalid tree extension."
    except ValueError:
        pass
