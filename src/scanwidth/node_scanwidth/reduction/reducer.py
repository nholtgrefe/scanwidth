"""Reducer abstraction for node-scanwidth algorithms."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
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

        block_infos: List[Tuple[Set, bool]]
        if self.config.use_sblocks:
            block_infos = sblocks(graph)
        else:
            block_infos = [(set(graph.nodes()), True)]

        # Blocks are independent and can later be processed in parallel.
        if self.config.parallel_sblocks and len(block_infos) > 1:
            block_results: List[Tuple[List, int]] = [([], 0)] * len(block_infos)
            with ThreadPoolExecutor(
                max_workers=self.config.sblock_max_workers,
            ) as executor:
                future_to_index = {
                    executor.submit(
                        self._solve_block,
                        graph,
                        block_info[0],
                        solver,
                        block_info[1],
                    ): i
                    for i, block_info in enumerate(block_infos)
                }
                for future in as_completed(future_to_index):
                    block_results[future_to_index[future]] = future.result()
        else:
            block_results = [
                self._solve_block(
                    graph, block_info[0], solver, block_info[1],
                )
                for block_info in block_infos
            ]

        sigma: List = []
        nsw = 0
        for partial_sigma, block_nsw in block_results:
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
        is_root_block: bool,
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

        # NOTE: Reticulation-path shortcut is intentionally disabled for now.
        # if self.config.use_reticulation_path_shortcut and not is_root_block:
        #     reticulation_path = self._solve_reticulation_path_block(subgraph)
        #     if reticulation_path is not None:
        #         partial_sigma, block_nsw = reticulation_path
        #         if self.config.use_chain_suppression:
        #             partial_sigma = self._unsuppress_chain_vertices(
        #                 partial_sigma, history,
        #             )
        #         return partial_sigma, block_nsw

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
        """Suppress chain parent in each chain-parent/chain-child pair.

        """
        reduced = subgraph.copy()
        history: List[Tuple[object, object]] = []

        changed = True
        while changed:
            changed = False
            for v in list(reduced.nodes()):
                if v not in reduced.nodes() or v not in chain_vertices(reduced):
                    continue
                w = list(reduced.successors(v))[0]
                if w not in chain_vertices(reduced):
                    continue
                u = list(reduced.predecessors(v))[0]
                reduced.remove_node(v)
                reduced.add_edge(u, w)
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

    # @staticmethod
    # def _solve_reticulation_path_block(
    #     subgraph: nx.DiGraph,
    # ) -> Optional[Tuple[List, int]]:
    #     """Solve special reticulation-path block using closed-form bound.

    #     Returns ``None`` when the structural preconditions are not met.
    #     """
        
    #     # By construction, this is biconnected.
    #     # if not nx.is_biconnected(subgraph.to_undirected()):
    #     #     return None

    #     reticulations = [v for v in subgraph.nodes() if subgraph.in_degree(v) >= 2]
    #     if not reticulations:
    #         return None

    #     topo_order = list(nx.topological_sort(subgraph))
    #     topo_pos = {v: i for i, v in enumerate(topo_order)}
    #     q = [v for v in topo_order if subgraph.in_degree(v) >= 2]
    #     for v in q:
    #         if subgraph.out_degree(v) > 1:
    #             return None
    #     for i in range(len(q) - 1):
    #         if (q[i], q[i + 1]) not in subgraph.edges():
    #             return None

    #     ancestor_sets: List[Set] = [set(nx.ancestors(subgraph, qi)) | {qi} for qi in q]
    #     A_sets: List[Set] = []
    #     seen: Set = set()
    #     for anc in ancestor_sets:
    #         a_i = anc.difference(seen)
    #         A_sets.append(a_i)
    #         seen = seen.union(anc)
    #     if seen != set(subgraph.nodes()):
    #         return None

    #     suffix_unions: List[Set] = [set() for _ in A_sets]
    #     running_suffix: Set = set()
    #     for i in reversed(range(len(A_sets))):
    #         running_suffix = running_suffix.union(A_sets[i])
    #         suffix_unions[i] = running_suffix.copy()

    #     Z_sets: List[Set] = []
    #     prefix_union: Set = set()
    #     for i, A_i in enumerate(A_sets):
    #         q_i = q[i]
    #         prefix_union = prefix_union.union(A_i)
    #         suffix_with_q_i = {q_i}.union(suffix_unions[i])
    #         parents_of_suffix = {
    #             u
    #             for (u, v) in subgraph.edges()
    #             if v in suffix_with_q_i and u not in suffix_with_q_i
    #         }
    #         Z_i = prefix_union.intersection(parents_of_suffix)
    #         Z_sets.append(Z_i)

    #     proof_order: List = []
    #     for A_i in A_sets:
    #         proof_order.extend(sorted(A_i, key=lambda x: topo_pos[x]))
    #     ordering = list(reversed(proof_order))

    #     extension = Extension(DAG(subgraph), ordering)
    #     nsw = max(len(Z_i) for Z_i in Z_sets)
    #     if nsw != extension.node_scanwidth():
    #         return None
    #     return ordering, nsw

