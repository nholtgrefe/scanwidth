"""Reducer abstraction for node-scanwidth algorithms."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Set, Tuple

import networkx as nx

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension
from scanwidth.node_scanwidth.reduction.config import ReducerConfig
from scanwidth.node_scanwidth.solver.base import Solver


@dataclass(frozen=True)
class Reducer:
    """Apply node-scanwidth reductions and delegate each block to a solver."""

    config: ReducerConfig = ReducerConfig()

    def reduce_and_solve(self, dag: DAG, solver: Solver) -> SolverResult:
        """Solve node scanwidth via reduction pipeline."""
        graph = dag.graph
        if self.config.use_tree_shortcut:
            tree = self._solve_tree(graph)
            if tree is not None:
                return SolverResult(value=tree[1], extension=Extension(dag, tree[0]))

        sblock_sets: List[Set]
        if self.config.use_sblocks:
            sblock_sets = self._sblocks(graph)
        else:
            sblock_sets = [set(graph.nodes())]

        # Blocks are independent and can later be processed in parallel.
        sigma: List = []
        nsw = 0
        for block_set in sblock_sets:
            partial_sigma, block_nsw = self._solve_block(graph, block_set, solver)
            nsw = max(nsw, block_nsw)
            partial_sigma = [v for v in partial_sigma if v not in sigma]
            sigma = partial_sigma + sigma

        extension = Extension(dag, sigma)
        return SolverResult(value=nsw, extension=extension)

    def _solve_block(
        self,
        graph: nx.DiGraph,
        sblock_set: Set,
        solver: Solver,
    ) -> Tuple[List, int]:
        """Solve one block and return partial order and block node-scanwidth."""
        if self.config.use_single_edge_shortcut:
            single_edge = self._solve_single_edge_block(graph, sblock_set)
            if single_edge is not None:
                return single_edge

        subgraph = graph.subgraph(sblock_set).copy()
        if self.config.use_tree_shortcut:
            tree = self._solve_tree(subgraph)
            if tree is not None:
                return tree

        history: List[Tuple[object, object]] = []
        if self.config.use_chain_suppression:
            subgraph, history = self._suppress_chain_vertices(subgraph)

        sub_result = solver.solve(DAG(subgraph))
        partial_sigma = list(sub_result.extension.ordering)
        if self.config.use_chain_suppression:
            partial_sigma = self._unsuppress_chain_vertices(partial_sigma, history)
        return partial_sigma, sub_result.value

    @staticmethod
    def _solve_single_edge_block(
        graph: nx.DiGraph,
        sblock_set: Set,
    ) -> Optional[Tuple[List, int]]:
        """Solve a block that is exactly one directed edge."""
        if len(sblock_set) != 2:
            return None
        u, v = sblock_set
        if (u, v) in graph.edges():
            return [v, u], 1
        return [u, v], 1

    @staticmethod
    def _solve_tree(graph: nx.DiGraph) -> Optional[Tuple[List, int]]:
        """Solve graph directly if it is a directed rooted tree."""
        if not Reducer._is_directed_tree(graph):
            return None
        sigma = list(reversed(list(nx.topological_sort(graph))))
        value = 0 if graph.number_of_nodes() <= 1 else 1
        return sigma, value

    @staticmethod
    def _is_chain_vertex(graph: nx.DiGraph, vertex: object) -> bool:
        """Return whether ``vertex`` is a chain vertex (in=1, out=1)."""
        return graph.in_degree(vertex) == 1 and graph.out_degree(vertex) == 1

    @staticmethod
    def _suppress_chain_vertices(
        subgraph: nx.DiGraph,
    ) -> Tuple[nx.DiGraph, List[Tuple[object, object]]]:
        """Suppress chain parent in each chain-parent/chain-child pair."""
        reduced = subgraph.copy()
        history: List[Tuple[object, object]] = []

        changed = True
        while changed:
            changed = False
            for v in list(reduced.nodes()):
                if v not in reduced.nodes() or not Reducer._is_chain_vertex(reduced, v):
                    continue
                w = list(reduced.successors(v))[0]
                if not Reducer._is_chain_vertex(reduced, w):
                    continue
                parent = list(reduced.predecessors(v))[0]
                reduced.remove_node(v)
                reduced.add_edge(parent, w)
                history.append((w, v))
                changed = True
                break

        return reduced, history

    @staticmethod
    def _unsuppress_chain_vertices(
        partial_sigma: List,
        history: List[Tuple[object, object]],
    ) -> List:
        """Undo chain-vertex suppression in reverse order."""
        restored = list(partial_sigma)
        for w, v in reversed(history):
            idx = restored.index(w)
            restored.insert(idx + 1, v)
        return restored

    @staticmethod
    def _is_directed_tree(graph: nx.DiGraph) -> bool:
        """Return whether graph is a directed rooted tree."""
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
        if len(roots) != 1:
            return False
        for v in graph.nodes():
            indeg = graph.in_degree(v)
            if v == roots[0]:
                if indeg != 0:
                    return False
            elif indeg != 1:
                return False
        return True

    @staticmethod
    def _sblocks(graph: nx.DiGraph) -> List[Set]:
        """Return node sets for s-blocks in merge order."""
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
