"""Reducer abstraction for edge-scanwidth algorithms."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import List, Optional, Set, Tuple

import networkx as nx

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.reduction.config import ReducerConfig
from scanwidth.edge_scanwidth.solver.base import Solver
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension
from scanwidth._utils import chain_vertices, sblocks


@dataclass(frozen=True)
class Reducer:
    """Apply s-block reduction and delegate each subproblem to a solver."""

    config: ReducerConfig = ReducerConfig()

    def reduce_and_solve(self, dag: DAG, solver: Solver) -> SolverResult:
        """Solve edge scanwidth via s-block decomposition.

        Splits ``dag`` into s-blocks, handles trivial blocks directly,
        contracts chain vertices inside non-trivial blocks, and delegates the
        contracted subproblem to ``solver``. Partial extensions are then
        reassembled into a single extension.

        Parameters
        ----------
        dag : DAG
            Input graph instance.
        solver : Solver
            Solver instance used on each contracted, non-trivial block.

        Returns
        -------
        SolverResult
            Solver result for the full graph.
        """
        graph = dag.graph
        if graph.number_of_nodes() == 1:
            (vertex,) = tuple(graph.nodes())
            return SolverResult(value=0, extension=Extension(dag, [vertex]))

        if self.config.use_sblocks:
            block_infos = sblocks(graph)
            sblock_sets = [block for block, _ in block_infos]
        else:
            sblock_sets = [set(graph.nodes())]
        if self.config.parallel_sblocks and len(sblock_sets) > 1:
            block_results: List[Tuple[List, int]] = [([], 0)] * len(sblock_sets)
            with ThreadPoolExecutor(
                max_workers=self.config.sblock_max_workers,
            ) as executor:
                future_to_index = {
                    executor.submit(self._solve_block, graph, sblock_set, solver): i
                    for i, sblock_set in enumerate(sblock_sets)
                }
                for future in as_completed(future_to_index):
                    block_results[future_to_index[future]] = future.result()
        else:
            block_results = [
                self._solve_block(graph, sblock_set, solver) for sblock_set in sblock_sets
            ]

        sigma: List = []
        sw = 0
        for partial_sigma, block_sw in block_results:
            sw = max(sw, block_sw)
            partial_sigma = [v for v in partial_sigma if v not in sigma]
            sigma = partial_sigma + sigma

        extension = Extension(dag, sigma)
        return SolverResult(value=sw, extension=extension)

    def _solve_block(
        self,
        graph: nx.DiGraph,
        sblock_set: Set,
        solver: Solver,
    ) -> tuple:
        """Solve a single s-block, returning (partial_sigma, block_sw)."""
        if self.config.use_single_edge_rule:
            single_edge = self._solve_single_edge_block(graph, sblock_set)
            if single_edge is not None:
                return single_edge

        subgraph = graph.subgraph(sblock_set)
        if self.config.use_single_root_cycle_rule:
            cycle = self._solve_single_root_cycle_block(subgraph)
            if cycle is not None:
                return cycle

        history: List[tuple] = []
        reduced_subgraph = subgraph.copy()
        if self.config.use_chain_suppression:
            reduced_subgraph, history = self._suppress_chain_vertices(subgraph)

        sub_result = solver.solve(DAG(reduced_subgraph))
        block_sw = sub_result.value
        partial_sigma = list(sub_result.extension.ordering)

        if self.config.use_chain_suppression:
            partial_sigma = self._unsuppress_chain_vertices(partial_sigma, history)

        return partial_sigma, block_sw

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
    def _solve_single_root_cycle_block(
        subgraph: nx.DiGraph,
    ) -> Optional[Tuple[List, int]]:
        """Solve a single-root cycle-like block with known scanwidth 2."""
        roots = [
            v for v in subgraph.nodes()
            if subgraph.in_degree(v) == 0 and subgraph.out_degree(v) == 2
        ]
        leafs = [
            v for v in subgraph.nodes()
            if subgraph.out_degree(v) == 0 and subgraph.in_degree(v) == 2
        ]
        chain_nodes = chain_vertices(subgraph)

        if (
            len(roots) == 1
            and len(leafs) == 1
            and len(chain_nodes) == len(subgraph.nodes()) - 2
        ):
            partial_sigma = list(reversed(list(nx.topological_sort(subgraph))))
            return partial_sigma, 2
        return None

    @staticmethod
    def _suppress_chain_vertices(
        subgraph: nx.DiGraph,
    ) -> Tuple[nx.DiGraph, List[tuple]]:
        """Contract suppressible chain vertices and return contraction history.

        Notes
        -----
        Chain vertices are the same as flow vertices (indegree 1, outdegree 1).
        """
        history: List[tuple] = []
        contracted = subgraph.copy()
        for v in chain_vertices(subgraph):
            u = list(contracted.predecessors(v))[0]
            w = list(contracted.successors(v))[0]
            if (u, w) not in contracted.edges():
                contracted.remove_node(v)
                contracted.add_edge(u, w)
                history.append((w, v))
        return contracted, history

    @staticmethod
    def _unsuppress_chain_vertices(partial_sigma: List, history: List[tuple]) -> List:
        """Undo suppressed chain-vertex contractions in reverse order."""
        restored = list(partial_sigma)
        for (w, v) in reversed(history):
            idx = restored.index(w)
            restored.insert(idx + 1, v)
        return restored

