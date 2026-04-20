"""Core correctness tests for scanwidth computations.

The tests in this module focus on small hand-crafted DAGs where expected
scanwidth values are easy to verify analytically.
"""

from __future__ import annotations

import networkx as nx

from scanwidth import DAG, Extension


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

    assert extension.is_extension()
    assert extension.scanwidth_at_vertex_i(0) == 2
    assert extension.scanwidth() == 2


def test_tree_extension_roundtrip_preserves_scanwidth() -> None:
    """Canonical tree extension preserves scanwidth and validity."""
    graph = _build_two_to_one()
    extension = Extension(graph=graph, sigma=[3, 1, 2])

    tree_extension = extension.canonical_tree_extension()
    roundtrip_extension = tree_extension.to_extension()

    assert tree_extension.is_canonical()
    assert tree_extension.scanwidth() == extension.scanwidth()
    assert roundtrip_extension.is_extension()
    assert roundtrip_extension.scanwidth() == extension.scanwidth()


def test_exact_methods_match_known_chain_scanwidth() -> None:
    """Exact methods compute scanwidth 1 on a directed chain."""
    dag = DAG(_build_chain())

    sw_exhaustive, ext_exhaustive = dag.optimal_scanwidth(reduced=False, method=1)
    sw_xp, ext_xp = dag.optimal_scanwidth(reduced=False, method=5)

    assert ext_exhaustive is not None
    assert ext_xp is not None
    assert sw_exhaustive == 1
    assert sw_xp == 1
    assert ext_exhaustive.scanwidth() == sw_exhaustive
    assert ext_xp.scanwidth() == sw_xp
    assert ext_exhaustive.is_extension()
    assert ext_xp.is_extension()


def test_heuristics_return_valid_extensions_and_consistent_scanwidth() -> None:
    """Heuristics return consistent scanwidth values and valid extensions."""
    dag = DAG(_build_chain())

    sw_greedy, ext_greedy = dag.greedy_heuristic(reduced=False)
    sw_cut, ext_cut = dag.cut_splitting_heuristic(reduced=False)
    sw_anneal, ext_anneal, _vals = dag.simulated_annealing(
        max_iter=5,
        verbose=False,
        reduced=False,
        seed=7,
    )

    assert sw_greedy == ext_greedy.scanwidth()
    assert sw_cut == ext_cut.scanwidth()
    assert sw_anneal == ext_anneal.scanwidth()
    assert ext_greedy.is_extension()
    assert ext_cut.is_extension()
    assert ext_anneal.is_extension()


def _assert_star_to_sink_scanwidth(expected_scanwidth: int) -> None:
    """Exact algorithm returns known scanwidths for simple star DAGs.

    Parameters
    ----------
    expected_scanwidth : int
        Expected optimal scanwidth of the generated DAG.
    """
    dag = DAG(_build_star_to_sink(expected_scanwidth))
    sw, ext = dag.optimal_scanwidth(reduced=False, method=1)

    assert ext is not None
    assert sw == expected_scanwidth
    assert ext.scanwidth() == expected_scanwidth
    assert ext.is_extension()


def test_star_to_sink_scanwidth_3() -> None:
    """Exact algorithm returns scanwidth 3 for the 3-star DAG."""
    _assert_star_to_sink_scanwidth(3)


def test_star_to_sink_scanwidth_4() -> None:
    """Exact algorithm returns scanwidth 4 for the 4-star DAG."""
    _assert_star_to_sink_scanwidth(4)
