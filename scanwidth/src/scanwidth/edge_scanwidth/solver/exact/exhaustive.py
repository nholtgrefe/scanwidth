"""Exhaustive-search exact solver for edge scanwidth."""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.solver.base import Solver
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension


@dataclass(frozen=True)
class ExhaustiveSolver(Solver):
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
        extensions = [ext[::-1] for ext in nx.all_topological_sorts(graph)]
        sw_values = [Extension(dag, sigma).edge_scanwidth() for sigma in extensions]
        best_sw = min(sw_values)
        best_sigma = extensions[sw_values.index(best_sw)]
        return SolverResult(value=best_sw, extension=Extension(dag, best_sigma))
