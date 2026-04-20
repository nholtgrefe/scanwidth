"""Exhaustive-search exact solver for node scanwidth."""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension
from scanwidth.node_scanwidth.solver.base import Solver


@dataclass(frozen=True)
class ExhaustiveSolver(Solver):
    """Exhaustive search over all topological orders for node scanwidth."""

    def solve(self, dag: DAG) -> SolverResult:
        """Solve node scanwidth by exhaustive search."""
        graph = dag.graph
        best_sw = None
        best_sigma = None
        for top_order in nx.all_topological_sorts(graph):
            sigma = list(reversed(list(top_order)))
            sw = Extension(dag, sigma).node_scanwidth()
            if best_sw is None or sw < best_sw:
                best_sw = sw
                best_sigma = sigma

        if best_sw is None or best_sigma is None:
            raise RuntimeError("No topological order found for DAG.")
        return SolverResult(value=best_sw, extension=Extension(dag, best_sigma))
