"""Random extension solver for edge scanwidth."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Set

import numpy as np

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.solver.base import Solver
from scanwidth._utils import delta_in
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension


@dataclass(frozen=True)
class RandomSolver(Solver):
    """Random-extension edge-scanwidth solver.

    Parameters
    ----------
    seed : int, optional
        Random seed for reproducibility.
    """

    seed: int = 42

    def solve(self, dag: DAG) -> SolverResult:
        """Generate a random extension and compute its edge scanwidth.

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
        rng = np.random.RandomState(self.seed)

        remaining: Set = set(graph.nodes())
        sigma: List = []
        sw = 0
        components_T: List[Set] = []

        while remaining:
            sub = graph.subgraph(remaining)
            leafs = [v for v in remaining if sub.out_degree(v) == 0]
            leaf = rng.choice(leafs)

            children = set(graph.successors(leaf))
            connected_vertices: Set = set()
            for comp in components_T:
                if not children.isdisjoint(comp):
                    connected_vertices = connected_vertices | comp
            connected_vertices.add(leaf)

            leaf_sw = delta_in(graph, connected_vertices)

            remaining.remove(leaf)
            sigma.append(leaf)
            sw = max(sw, leaf_sw)

            components_T = [
                comp for comp in components_T
                if not comp.issubset(connected_vertices)
            ]
            components_T.append(connected_vertices)

        extension = Extension(dag, sigma)
        return SolverResult(value=sw, extension=extension)
