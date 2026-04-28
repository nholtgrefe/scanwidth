"""Core correctness tests for scanwidth computations.

The tests in this module focus on small hand-crafted DAGs where expected
scanwidth values are easy to verify analytically.
"""

from __future__ import annotations

import builtins
import networkx as nx
import pytest  # pyright: ignore[reportMissingImports]

from scanwidth import DAG, Extension, TreeExtension, edge_scanwidth, node_scanwidth
from scanwidth.edge_scanwidth.reduction.config import ReducerConfig as EdgeReducerConfig
from scanwidth.node_scanwidth.reduction.config import ReducerConfig as NodeReducerConfig
from scanwidth.node_scanwidth.reduction.reducer import Reducer
from scanwidth._utils import (
    chain_vertices,
    delta_in,
    find_component,
    infinity_for,
    is_directed_tree,
    sblocks,
)


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
        dag, algorithm="brute_force", reduce=False,
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


def test_node_scanwidth_algorithms_on_chain() -> None:
    """Node-scanwidth algorithms return width 1 on a directed chain."""
    dag = DAG(_build_chain())

    sw_xp, ext_xp = node_scanwidth(dag, algorithm="xp", reduce=False)
    sw_exhaustive, ext_exhaustive = node_scanwidth(dag, algorithm="brute_force")
    sw_greedy, ext_greedy = node_scanwidth(dag, algorithm="greedy")
    sw_random, ext_random = node_scanwidth(dag, algorithm="random", seed=7)
    sw_anneal, ext_anneal = node_scanwidth(
        dag,
        algorithm="simulated_annealing",
        max_iter=5,
        verbose=False,
        seed=7,
    )

    assert sw_xp == 1
    assert sw_exhaustive == 1
    assert sw_greedy == 1
    assert sw_random == 1
    assert sw_anneal == 1
    assert ext_xp.node_scanwidth() == 1
    assert ext_exhaustive.node_scanwidth() == 1
    assert ext_greedy.node_scanwidth() == 1
    assert ext_random.node_scanwidth() == 1
    assert ext_anneal.node_scanwidth() == 1


def test_node_scanwidth_exhaustive_matches_known_two_to_one_value() -> None:
    """Exhaustive node-scanwidth solver matches known optimum."""
    dag = DAG(_build_two_to_one())
    sw, extension = node_scanwidth(dag, algorithm="brute_force")
    assert sw == 2
    assert extension.node_scanwidth() == 2


def test_node_scanwidth_reducer_tree_rule_returns_one() -> None:
    """Node reducer applies tree shortcut with node scanwidth 1."""
    graph = nx.DiGraph([(1, 2), (1, 3), (2, 4), (2, 5)])
    dag = DAG(graph)
    value, extension = node_scanwidth(dag, algorithm="greedy", reduce=True)
    assert value == 1
    assert extension.node_scanwidth() == 1


def test_node_scanwidth_singleton_graph_shortcut_returns_zero() -> None:
    """Node reducer returns zero and trivial extension on singleton graph."""
    graph = nx.DiGraph()
    graph.add_node("only")
    dag = DAG(graph)
    value, extension = node_scanwidth(dag, algorithm="greedy", reduce=True)
    assert value == 0
    assert extension.ordering == ["only"]
    assert extension.node_scanwidth() == 0


def test_edge_scanwidth_singleton_graph_shortcut_returns_zero() -> None:
    """Edge reducer returns zero and trivial extension on singleton graph."""
    graph = nx.DiGraph()
    graph.add_node("only")
    dag = DAG(graph)
    value, extension = edge_scanwidth(dag, algorithm="greedy", reduce=True)
    assert value == 0
    assert extension.ordering == ["only"]
    assert extension.edge_scanwidth() == 0


def test_edge_scanwidth_parallel_sblocks_matches_sequential() -> None:
    """Parallel and sequential edge s-block solving return same value."""
    dag = DAG(_build_chain())
    seq_value, seq_ext = edge_scanwidth(
        dag,
        algorithm="brute_force",
        reduce=True,
        reducer_config=EdgeReducerConfig(parallel_sblocks=False),
    )
    par_value, par_ext = edge_scanwidth(
        dag,
        algorithm="brute_force",
        reduce=True,
        reducer_config=EdgeReducerConfig(parallel_sblocks=True, sblock_max_workers=2),
    )
    assert par_value == seq_value
    assert par_value == par_ext.edge_scanwidth() == seq_ext.edge_scanwidth()


def test_node_scanwidth_parallel_sblocks_matches_sequential() -> None:
    """Parallel and sequential node s-block solving return same value."""
    dag = DAG(_build_chain())
    seq_value, seq_ext = node_scanwidth(
        dag,
        algorithm="brute_force",
        reduce=True,
        reducer_config=NodeReducerConfig(parallel_sblocks=False),
    )
    par_value, par_ext = node_scanwidth(
        dag,
        algorithm="brute_force",
        reduce=True,
        reducer_config=NodeReducerConfig(parallel_sblocks=True, sblock_max_workers=2),
    )
    assert par_value == seq_value
    assert par_value == par_ext.node_scanwidth() == seq_ext.node_scanwidth()


def test_node_scanwidth_reticulation_path_shortcut_matches_bruteforce() -> None:
    """Reticulation-path shortcut keeps brute-force optimum on diamond DAG."""
    graph = nx.DiGraph([("r", "a"), ("r", "b"), ("a", "q"), ("b", "q")])
    dag = DAG(graph)

    value_reduced, ext_reduced = node_scanwidth(
        dag, algorithm="brute_force", reduce=True,
    )
    value_direct, ext_direct = node_scanwidth(
        dag, algorithm="brute_force", reduce=False,
    )
    assert value_reduced == value_direct == 2
    assert ext_reduced.node_scanwidth() == ext_direct.node_scanwidth() == 2


def test_node_scanwidth_reticulation_path_shortcut_toggle() -> None:
    """Reticulation-path shortcut toggle keeps optimal node scanwidth value."""
    graph = nx.DiGraph([("r", "a"), ("r", "b"), ("a", "q"), ("b", "q")])
    dag = DAG(graph)

    value_on, ext_on = node_scanwidth(
        dag,
        algorithm="brute_force",
        reduce=True,
        reducer_config=NodeReducerConfig(use_reticulation_path_shortcut=True),
    )
    value_off, ext_off = node_scanwidth(
        dag,
        algorithm="brute_force",
        reduce=True,
        reducer_config=NodeReducerConfig(use_reticulation_path_shortcut=False),
    )
    assert value_on == value_off == 2
    assert ext_on.node_scanwidth() == ext_off.node_scanwidth() == 2


def test_node_scanwidth_reducer_chain_suppression_rule() -> None:
    """Reducer suppresses chain-parent in a chain-parent/chain-child pair."""
    # 0 -> 1 -> 2 -> 3, where 1 and 2 are chain vertices and 1 is parent of 2.
    subgraph = nx.DiGraph([(0, 1), (1, 2), (2, 3)])
    reduced, history = Reducer._suppress_chain_vertices(subgraph)
    assert set(reduced.nodes()) == {0, 2, 3}
    assert (0, 2) in reduced.edges()
    restored = Reducer._unsuppress_chain_vertices([3, 2, 0], history)
    assert restored == [3, 2, 1, 0]


def test_global_sblocks_returns_expected_blocks_for_chain() -> None:
    """Global s-block helper returns expected biconnected blocks on a chain."""
    graph = _build_chain()
    block_infos = sblocks(graph)
    blocks = [block for block, _ in block_infos]
    is_root_block = [flag for _, flag in block_infos]
    assert blocks == [{1, 2}, {2, 3}]
    assert is_root_block == [False, False]


def test_global_sblocks_marks_special_root_block_for_multi_root_graph() -> None:
    """S-block helper marks root block when graph has multiple roots."""
    graph = nx.DiGraph([(1, 3), (2, 3)])
    block_infos = sblocks(graph)
    blocks = [block for block, _ in block_infos]
    is_root_block = [flag for _, flag in block_infos]
    assert len(blocks) == 1
    assert blocks[0] == {1, 2, 3}
    assert is_root_block == [True]


def test_global_edge_solver_utils_match_basic_expectations() -> None:
    """Global edge helper utilities keep expected semantics."""
    graph = _build_two_to_one()
    vertices = {3}
    components = [{1, 2}, {3}]

    assert delta_in(graph, vertices) == 2
    assert find_component(components, 3) == {3}
    assert find_component(components, "missing") == set()
    assert infinity_for(graph) == graph.number_of_edges() + 1


def test_global_chain_vertices_matches_flow_vertex_definition() -> None:
    """Global chain-vertex helper returns degree-(1,1) vertices."""
    graph = nx.DiGraph([(0, 1), (1, 2), (2, 3), (2, 4)])
    assert chain_vertices(graph) == [1]


def test_global_is_directed_tree_detects_tree_shape() -> None:
    """Global directed-tree helper identifies directed rooted trees."""
    tree_graph = nx.DiGraph([(1, 2), (1, 3), (3, 4)])
    non_tree_graph = nx.DiGraph([(1, 2), (3, 2)])
    assert is_directed_tree(tree_graph)
    assert not is_directed_tree(non_tree_graph)


def test_node_scanwidth_ilp_matches_exhaustive_on_chain() -> None:
    """ILP node-scanwidth solver matches exhaustive optimum on a chain."""
    dag = DAG(_build_chain())
    sw_ilp, ext_ilp = node_scanwidth(dag, algorithm="ilp")
    sw_ilp_backend_scipy, ext_ilp_backend_scipy = node_scanwidth(
        dag,
        algorithm="ilp",
        backend="scipy",
    )
    sw_exact, ext_exact = node_scanwidth(dag, algorithm="brute_force")
    assert sw_ilp == sw_exact == 1
    assert sw_ilp_backend_scipy == sw_exact == 1
    assert ext_ilp.node_scanwidth() == ext_exact.node_scanwidth() == 1
    assert ext_ilp_backend_scipy.node_scanwidth() == ext_exact.node_scanwidth() == 1


def test_node_scanwidth_ilp_matches_exhaustive_on_two_to_one() -> None:
    """ILP node-scanwidth solver matches exhaustive optimum on two-to-one DAG."""
    dag = DAG(_build_two_to_one())
    sw_ilp, ext_ilp = node_scanwidth(dag, algorithm="ilp")
    sw_ilp_backend_scipy, ext_ilp_backend_scipy = node_scanwidth(
        dag,
        algorithm="ilp",
        backend="scipy",
    )
    sw_exact, ext_exact = node_scanwidth(dag, algorithm="brute_force")
    assert sw_ilp == sw_exact == 2
    assert sw_ilp_backend_scipy == sw_exact == 2
    assert ext_ilp.node_scanwidth() == ext_exact.node_scanwidth() == 2
    assert ext_ilp_backend_scipy.node_scanwidth() == ext_exact.node_scanwidth() == 2


def test_node_scanwidth_ilp_gurobi_matches_exhaustive_on_chain() -> None:
    """Gurobi-backend ILP solver matches exhaustive optimum on a chain."""
    pytest.importorskip("gurobipy")
    dag = DAG(_build_chain())
    sw_ilp, ext_ilp = node_scanwidth(dag, algorithm="ilp", backend="gurobi")
    sw_exact, ext_exact = node_scanwidth(dag, algorithm="brute_force")
    assert sw_ilp == sw_exact == 1
    assert ext_ilp.node_scanwidth() == ext_exact.node_scanwidth() == 1


def test_node_scanwidth_ilp_rejects_unknown_backend() -> None:
    """ILP node API rejects unsupported backend names."""
    dag = DAG(_build_chain())
    with pytest.raises(ValueError, match="Unknown ILP backend"):
        _ = node_scanwidth(dag, algorithm="ilp", backend="unknown", reduce=False)


def test_node_scanwidth_ilp_scipy_backend_missing_dependency_raises_import_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """SciPy backend raises informative ImportError when SciPy is unavailable."""
    original_import = builtins.__import__

    def _patched_import(name: str, *args: object, **kwargs: object) -> object:
        if name.startswith("scipy"):
            raise ImportError("mocked missing scipy")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _patched_import)
    dag = DAG(_build_chain())
    with pytest.raises(ImportError, match="scanwidth\\[scipy\\]"):
        _ = node_scanwidth(dag, algorithm="ilp", backend="scipy", reduce=False)


def test_node_scanwidth_ilp_gurobi_backend_missing_dependency_raises_import_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Gurobi backend raises informative ImportError when gurobipy is unavailable."""
    original_import = builtins.__import__

    def _patched_import(name: str, *args: object, **kwargs: object) -> object:
        if name == "gurobipy":
            raise ImportError("mocked missing gurobipy")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _patched_import)
    dag = DAG(_build_chain())
    with pytest.raises(ImportError, match="scanwidth\\[gurobi\\]"):
        _ = node_scanwidth(dag, algorithm="ilp", backend="gurobi", reduce=False)


def test_node_scanwidth_xp_matches_exhaustive_on_two_to_one() -> None:
    """XP node-scanwidth solver matches exhaustive optimum on two-to-one DAG."""
    dag = DAG(_build_two_to_one())
    sw_xp, ext_xp = node_scanwidth(dag, algorithm="xp", reduce=False)
    sw_exact, ext_exact = node_scanwidth(dag, algorithm="brute_force", reduce=False)
    assert sw_xp == sw_exact == 2
    assert ext_xp.node_scanwidth() == ext_exact.node_scanwidth() == 2


def test_node_scanwidth_entrypoint_xp_rejects_fixed_k() -> None:
    """Public XP node entrypoint rejects fixed-k mode."""
    dag = DAG(_build_chain())
    try:
        _ = node_scanwidth(dag, algorithm="xp", k=1, reduce=False)
        assert False, "Expected ValueError when passing k to public node xp API."
    except ValueError:
        pass


def _assert_star_to_sink_scanwidth(expected_scanwidth: int) -> None:
    """Exact algorithm returns known scanwidths for simple star DAGs.

    Parameters
    ----------
    expected_scanwidth : int
        Expected optimal scanwidth of the generated DAG.
    """
    dag = DAG(_build_star_to_sink(expected_scanwidth))
    sw, ext = edge_scanwidth(dag, algorithm="brute_force", reduce=False)

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


