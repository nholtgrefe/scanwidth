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
from scanwidth._utils import chain_vertices, is_directed_tree, sblocks


@dataclass(frozen=True)
class Reducer:
    """Apply node-scanwidth reductions and delegate each block to a solver."""

    config: ReducerConfig = ReducerConfig()

    def reduce_and_solve(self, dag: DAG, solver: Solver) -> SolverResult:
        """Solve node scanwidth via reduction pipeline."""
        graph = dag.graph
        if graph.number_of_nodes() == 1:
            (vertex,) = tuple(graph.nodes())
            return SolverResult(value=0, extension=Extension(dag, [vertex]))

        if self.config.use_tree_shortcut:
            tree = self._solve_tree(graph)
            if tree is not None:
                return SolverResult(value=tree[1], extension=Extension(dag, tree[0]))

        sblock_sets: List[Set]
        if self.config.use_sblocks:
            sblock_sets = sblocks(graph)
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
        if not is_directed_tree(graph):
            return None
        sigma = list(reversed(list(nx.topological_sort(graph))))
        value = 0 if graph.number_of_nodes() <= 1 else 1
        return sigma, value

    @staticmethod
    def _suppress_chain_vertices(
        subgraph: nx.DiGraph,
    ) -> Tuple[nx.DiGraph, List[Tuple[object, object]]]:
        """Contract suppressible chain vertices and return contraction history.

        Notes
        -----
        Chain vertices are the same as flow vertices (indegree 1, outdegree 1).
        """
        reduced = subgraph.copy()
        history: List[Tuple[object, object]] = []

        for v in chain_vertices(subgraph):
            u = list(reduced.predecessors(v))[0]
            w = list(reduced.successors(v))[0]
            if (u, w) not in reduced.edges():
                reduced.remove_node(v)
                reduced.add_edge(u, w)
                history.append((w, v))

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

