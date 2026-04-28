"""Greedy solver for node scanwidth."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Set

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension
from scanwidth.node_scanwidth.solver.base import Solver
from scanwidth._utils import infinity_for
from scanwidth.node_scanwidth.solver._utils import node_bag_size


@dataclass(frozen=True)
class GreedySolver(Solver):
    """Greedy node-scanwidth solver."""

    def solve(self, dag: DAG) -> SolverResult:
        """Solve node scanwidth using the greedy heuristic."""
        graph = dag.graph
        infinity = infinity_for(graph)

        remaining: Set = set(graph.nodes())
        sigma: List = []
        sw = 0
        components_t: List[Set] = []

        while remaining:
            chosen = None
            chosen_sw = infinity
            sub = graph.subgraph(remaining)
            leafs = [v for v in remaining if sub.out_degree(v) == 0]

            connection_dict: Dict = {}
            for leaf in leafs:
                children = set(graph.successors(leaf))
                connected_vertices: Set = set()
                for comp in components_t:
                    if not children.isdisjoint(comp):
                        connected_vertices = connected_vertices | comp
                connected_vertices.add(leaf)
                connection_dict[leaf] = connected_vertices

                candidate_sw = node_bag_size(graph, connected_vertices)
                if candidate_sw < chosen_sw:
                    chosen = leaf
                    chosen_sw = candidate_sw

            remaining.remove(chosen)
            sigma.append(chosen)
            sw = max(sw, chosen_sw)

            components_t = [
                comp for comp in components_t if not comp.issubset(connection_dict[chosen])
            ]
            components_t.append(connection_dict[chosen])

        extension = Extension(dag, sigma)
        return SolverResult(value=sw, extension=extension)
