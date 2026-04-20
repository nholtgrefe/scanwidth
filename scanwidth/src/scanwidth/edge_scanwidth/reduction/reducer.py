"""Reducer abstraction for edge-scanwidth algorithms."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Set

import networkx as nx

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.solver.base import Solver
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension


@dataclass(frozen=True)
class Reducer:
    """Apply s-block reduction and delegate each subproblem to a solver."""

    def reduce_and_solve(self, dag: DAG, solver: Solver) -> SolverResult:
        """Solve edge scanwidth via s-block decomposition.

        Splits ``dag`` into s-blocks, handles trivial blocks directly,
        contracts flow-edges inside non-trivial blocks, and delegates the
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
        sblock_sets = self._sblocks(graph)
        sigma: List = []
        sw = 0

        for sblock_set in sblock_sets:
            partial_sigma, block_sw = self._solve_block(
                graph, sblock_set, solver,
            )
            sw = max(sw, block_sw)
            partial_sigma = [v for v in partial_sigma if v not in sigma]
            sigma = partial_sigma + sigma

        extension = Extension(dag, sigma)
        return SolverResult(value=sw, extension=extension)

    @staticmethod
    def _solve_block(
        graph: nx.DiGraph,
        sblock_set: Set,
        solver: Solver,
    ) -> tuple:
        """Solve a single s-block, returning (partial_sigma, block_sw)."""
        if len(sblock_set) == 2:
            u, v = sblock_set
            if (u, v) in graph.edges():
                return [v, u], 1
            return [u, v], 1

        subgraph = graph.subgraph(sblock_set)
        roots = [
            v for v in subgraph.nodes()
            if subgraph.in_degree(v) == 0 and subgraph.out_degree(v) == 2
        ]
        leafs = [
            v for v in subgraph.nodes()
            if subgraph.out_degree(v) == 0 and subgraph.in_degree(v) == 2
        ]
        flow_nodes = [
            v for v in subgraph.nodes()
            if subgraph.out_degree(v) == 1 and subgraph.in_degree(v) == 1
        ]

        if (
            len(roots) == 1
            and len(leafs) == 1
            and len(flow_nodes) == len(subgraph.nodes()) - 2
        ):
            partial_sigma = list(reversed(list(nx.topological_sort(subgraph))))
            return partial_sigma, 2

        history: List[tuple] = []
        contracted = subgraph.copy()
        for v in flow_nodes:
            u = list(contracted.predecessors(v))[0]
            w = list(contracted.successors(v))[0]
            if (u, w) not in contracted.edges():
                contracted.remove_node(v)
                contracted.add_edge(u, w)
                history.append((w, v))

        sub_result = solver.solve(DAG(contracted))
        block_sw = sub_result.value
        partial_sigma = list(sub_result.extension.ordering)

        history.reverse()
        for (w, v) in history:
            idx = partial_sigma.index(w)
            partial_sigma.insert(idx + 1, v)

        return partial_sigma, block_sw

    @staticmethod
    def _sblocks(graph: nx.DiGraph) -> List[Set]:
        """Return node sets for s-blocks in merge order.

        Uses a reversed DFS in the s-block-cut-tree starting at the
        rootblock, so the rootblock comes last.

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
                    (root1, root2) not in aux.edges()
                    and (root2, root1) not in aux.edges
                    and root1 != root2
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
