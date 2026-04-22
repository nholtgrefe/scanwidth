"""Random extension solver for node scanwidth."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Set

import numpy as np

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension
from scanwidth.node_scanwidth.solver.base import Solver
from scanwidth.node_scanwidth.solver._utils import node_bag_size


@dataclass(frozen=True)
class RandomSolver(Solver):
    """Random-extension node-scanwidth solver."""

    seed: int = 42

    def solve(self, dag: DAG) -> SolverResult:
        """Generate a random extension and compute its node scanwidth."""
        graph = dag.graph
        rng = np.random.RandomState(self.seed)

        remaining: Set = set(graph.nodes())
        sigma: List = []
        sw = 0
        components_t: List[Set] = []

        while remaining:
            sub = graph.subgraph(remaining)
            leafs = [v for v in remaining if sub.out_degree(v) == 0]
            leaf = rng.choice(leafs)

            children = set(graph.successors(leaf))
            connected_vertices: Set = set()
            for comp in components_t:
                if not children.isdisjoint(comp):
                    connected_vertices = connected_vertices | comp
            connected_vertices.add(leaf)

            leaf_sw = node_bag_size(graph, connected_vertices)

            remaining.remove(leaf)
            sigma.append(leaf)
            sw = max(sw, leaf_sw)

            components_t = [
                comp for comp in components_t if not comp.issubset(connected_vertices)
            ]
            components_t.append(connected_vertices)

        extension = Extension(dag, sigma)
        return SolverResult(value=sw, extension=extension)
