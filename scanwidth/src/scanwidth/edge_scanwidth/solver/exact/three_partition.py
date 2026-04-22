"""3-partition recursion exact solver for edge scanwidth."""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import List, Set, Tuple

import networkx as nx

from scanwidth.dag import DAG
from scanwidth.edge_scanwidth.solver.base import Solver
from scanwidth._utils import (
    delta_in,
    find_component,
    infinity_for,
)
from scanwidth.edge_scanwidth.types import SolverResult
from scanwidth.extension import Extension


@dataclass(frozen=True)
class ThreePartitionSolver(Solver):
    """Exact solver based on ordered 3-partition recursion.

    Recursively splits the vertex set into an ordered 3-partition
    ``(L, W, R)`` with ``W`` halved at each step, under the constraint
    that there are no edges from ``W'`` to ``W \\ W'``.
    """

    def solve(self, dag: DAG) -> SolverResult:
        """Solve edge scanwidth using 3-partition recursion.

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
        infinity = infinity_for(graph)
        sw, sigma = self._partial_scanwidth(
            graph, set(), set(graph.nodes()), set(), infinity,
        )
        return SolverResult(value=sw, extension=Extension(dag, sigma))

    def _partial_scanwidth(
        self,
        graph: nx.DiGraph,
        L: Set,
        W: Set,
        R: Set,
        infinity: int,
    ) -> Tuple[int, List]:
        """Compute partial scanwidth for the ordered partition ``(L, W, R)``."""
        psw = infinity
        sigma: List = []

        if len(W) == 1:
            (w,) = W
            subgraph = graph.subgraph(W.union(L))
            components = list(nx.weakly_connected_components(subgraph))
            component = find_component(components, w)
            psw = delta_in(graph, component)
            sigma = [w]
        elif len(W) > 1:
            size = len(W) // 2
            for W_prime in map(set, itertools.combinations(W, size)):
                check = len([
                    (u, v) for u in W_prime
                    for v in W.difference(W_prime)
                    if (u, v) in graph.edges()
                ])
                if check != 0:
                    continue

                psw_1, sigma_1 = self._partial_scanwidth(
                    graph, L, W_prime,
                    R.union(W.difference(W_prime)),
                    infinity,
                )
                psw_2, sigma_2 = self._partial_scanwidth(
                    graph,
                    L.union(W_prime),
                    W.difference(W_prime),
                    R,
                    infinity,
                )
                psw_prime = max(psw_1, psw_2)
                if psw_prime < psw:
                    psw = psw_prime
                    sigma = sigma_1 + sigma_2

        return psw, sigma
