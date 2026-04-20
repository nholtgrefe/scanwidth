"""Cut-splitting solver for edge scanwidth."""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import List, Optional, Set, Tuple

import networkx as nx
from networkx.algorithms.flow import shortest_augmenting_path

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.solver.base import Solver
from scanwidth.edge_scanwidth.solver.utils import infinity_for
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension


@dataclass(frozen=True)
class CutSplittingSolver(Solver):
    """Recursive cut-splitting edge-scanwidth solver."""

    def solve(self, dag: DAG) -> SolverResult:
        """Solve edge scanwidth using recursive cut-splitting.

        Parameters
        ----------
        dag : DAG
            Input graph instance.

        Returns
        -------
        SolverResult
            Standardized solver result.
        """
        graph = dag.graph
        sigma = self._recursive_cut_splitting(graph)
        extension = Extension(dag, sigma)
        return SolverResult(value=extension.edge_scanwidth(), extension=extension)

    def _recursive_cut_splitting(self, graph: nx.DiGraph) -> List:
        """Recursively split the graph along the minimum DAG cut."""
        edge_cut, source_set, sink_set = self._minimum_DAG_cut(
            graph, non_trivial=True,
        )

        if edge_cut is None:
            return list(reversed(list(nx.topological_sort(graph))))

        subgraph_top = graph.subgraph(source_set).copy()
        subgraph_down = graph.subgraph(sink_set).copy()

        super_leaf = "super_leaf_" + str(random.getrandbits(64))
        super_root = "super_root_" + str(random.getrandbits(64))
        subgraph_top.add_node(super_leaf, merged=True)
        subgraph_down.add_node(super_root, merged=True)

        for (u, v) in edge_cut:
            if "weight" not in graph[u][v]:
                weight = 1
            else:
                weight = graph[u][v]["weight"]

            if subgraph_top.has_edge(u, super_leaf):
                subgraph_top[u][super_leaf]["weight"] += weight
            else:
                subgraph_top.add_edge(u, super_leaf, weight=weight)

            if subgraph_down.has_edge(super_root, v):
                subgraph_down[super_root][v]["weight"] += weight
            else:
                subgraph_down.add_edge(super_root, v, weight=weight)

        sigma_top = self._recursive_cut_splitting(subgraph_top)
        sigma_down = self._recursive_cut_splitting(subgraph_down)
        sigma_top.remove(super_leaf)
        sigma_down.remove(super_root)

        return sigma_down + sigma_top

    @staticmethod
    def _minimum_DAG_cut(
        graph: nx.DiGraph,
        non_trivial: bool = True,
    ) -> Tuple[Optional[List], Optional[Set], Optional[Set]]:
        """Find the smallest DAG cut of ``graph``.

        Parameters
        ----------
        graph : nx.DiGraph
            Input graph.
        non_trivial : bool, optional
            If True, ignore cuts that only isolate a single merged node.

        Returns
        -------
        Tuple[Optional[List], Optional[Set], Optional[Set]]
            Edges in the cut, and the node sets on each side. Returns
            ``(None, None, None)`` when no cut exists.
        """
        infinity = infinity_for(graph)
        aux_graph = graph.copy()

        for (u, v) in aux_graph.edges():
            if "weight" not in aux_graph[u][v]:
                aux_graph[u][v]["weight"] = 1

        aux_graph.add_weighted_edges_from(
            [(v, u, infinity) for (u, v) in aux_graph.edges]
        )

        roots = [v for v in graph.nodes() if graph.in_degree(v) == 0]
        leafs = [v for v in graph.nodes() if graph.out_degree(v) == 0]

        cut = infinity
        source_set: Set = set()
        sink_set: Set = set()
        size_diff = infinity + 1

        if not non_trivial:
            for root in roots:
                for leaf in leafs:
                    a, (b, c) = nx.minimum_cut(
                        aux_graph,
                        root,
                        leaf,
                        capacity="weight",
                        flow_func=shortest_augmenting_path,
                    )
                    d = abs(len(b) - len(c))
                    if a < cut or (a == cut and d < size_diff):
                        cut, source_set, sink_set, size_diff = a, b, c, d
        else:
            for root in roots:
                if "merged" in graph.nodes[root]:
                    root_childs = list(graph.successors(root))
                else:
                    root_childs = [root]

                for leaf in leafs:
                    if "merged" in graph.nodes[leaf]:
                        leaf_parents = list(graph.predecessors(leaf))
                    else:
                        leaf_parents = [leaf]

                    for root_child in root_childs:
                        for leaf_parent in leaf_parents:
                            if root_child == leaf_parent:
                                continue
                            a, (b, c) = nx.minimum_cut(
                                aux_graph,
                                root_child,
                                leaf_parent,
                                capacity="weight",
                                flow_func=shortest_augmenting_path,
                            )
                            d = abs(len(b) - len(c))
                            if a < cut or (a == cut and d < size_diff):
                                cut, source_set, sink_set, size_diff = a, b, c, d

        if cut == infinity:
            return None, None, None

        edge_cut: List = []
        for (u, v) in graph.edges():
            if u in source_set and v in sink_set:
                edge_cut.append((u, v))

        return edge_cut, source_set, sink_set
