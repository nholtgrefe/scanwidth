"""Exhaustive-search exact solver for edge scanwidth."""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.solver.base import Solver
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension


@dataclass(frozen=True)
class BruteForceSolver(Solver):
    """Exhaustive search over all topological orders.

    Enumerates every topological order of the input DAG, evaluates the
    scanwidth of the corresponding extension, and returns the best one.
    Only tractable on very small inputs.
    """

    def solve(self, dag: DAG) -> SolverResult:
        """Solve edge scanwidth by exhaustive search.

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
        best_sw = None
        best_sigma = None
        for top_order in nx.all_topological_sorts(graph):
            sigma = list(reversed(list(top_order)))
            sw = Extension(dag, sigma).edge_scanwidth()
            if best_sw is None or sw < best_sw:
                best_sw = sw
                best_sigma = sigma

        if best_sw is None or best_sigma is None:
            raise RuntimeError("No topological order found for DAG.")
        return SolverResult(value=best_sw, extension=Extension(dag, best_sigma))
